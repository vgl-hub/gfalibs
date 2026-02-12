// ------------------------------
// Lightweight semaphore (hybrid)
// ------------------------------
class LightweightSemaphore {
public:
	explicit LightweightSemaphore(int initial = 0)
		: count_(initial), wakeups_(0) {}

	void acquire() {
		int old = count_.fetch_sub(1, std::memory_order_acq_rel);
		if (old > 0) return; // fast path

		std::unique_lock<std::mutex> lk(m_);
		cv_.wait(lk, [&]{ return wakeups_ > 0; });
		--wakeups_; // consume exactly one wakeup
	}

	bool try_acquire() {
		int c = count_.load(std::memory_order_relaxed);
		while (c > 0) {
			if (count_.compare_exchange_weak(
					c, c - 1,
					std::memory_order_acq_rel,
					std::memory_order_relaxed)) return true;
			// c reloaded on failure
		}
		return false;
	}

	void release() {
		int old = count_.fetch_add(1, std::memory_order_release);
		if (old < 0) { // sleepers present
			std::lock_guard<std::mutex> lk(m_);
			++wakeups_;
			cv_.notify_one();
		}
	}

private:
	std::atomic<int> count_;   // can go negative: -sleepers
	std::mutex m_;
	std::condition_variable cv_;
	int wakeups_;
};

// ---------------------------------------
// Lock-free bounded MPMC ring (Vyukov)
// ---------------------------------------
template<class T>
class MpmcRing {
public:
	explicit MpmcRing(size_t capacity) {
		n_ = 1;
		while (n_ < capacity) n_ <<= 1;
		mask_ = n_ - 1;
		slots_.resize(n_);
		for (size_t i = 0; i < n_; ++i)
			slots_[i].seq.store(i, std::memory_order_relaxed);
		head_.store(0, std::memory_order_relaxed);
		tail_.store(0, std::memory_order_relaxed);
	}

	// NOTE: rvalue ref here
	bool try_enqueue(T&& v) {
		size_t pos = tail_.load(std::memory_order_relaxed);
		for (;;) {
			Slot& s = slots_[pos & mask_];
			size_t seq = s.seq.load(std::memory_order_acquire);
			std::intptr_t diff =
				static_cast<std::intptr_t>(seq) -
				static_cast<std::intptr_t>(pos);

			if (diff == 0) {
				// slot is ready for us
				if (tail_.compare_exchange_weak(
						pos, pos + 1,
						std::memory_order_relaxed,
						std::memory_order_relaxed)) {
					// OWNERSHIP ACQUIRED → now it's safe to move
					s.val = std::move(v);
					s.seq.store(pos + 1, std::memory_order_release);
					return true;
				}
				// CAS failed => someone else claimed it; retry with new pos
			} else if (diff < 0) {
				// queue is full
				return false;
			} else {
				// another producer moved tail forward; reload
				pos = tail_.load(std::memory_order_relaxed);
			}
		}
	}

	bool try_dequeue(T& out) {
		size_t pos = head_.load(std::memory_order_relaxed);
		for (;;) {
			Slot& s = slots_[pos & mask_];
			size_t seq = s.seq.load(std::memory_order_acquire);
			std::intptr_t diff =
				static_cast<std::intptr_t>(seq) -
				static_cast<std::intptr_t>(pos + 1);

			if (diff == 0) {
				if (head_.compare_exchange_weak(
						pos, pos + 1,
						std::memory_order_relaxed,
						std::memory_order_relaxed)) {
					out = std::move(s.val);
					s.seq.store(pos + n_, std::memory_order_release);
					return true;
				}
			} else if (diff < 0) {
				// empty
				return false;
			} else {
				pos = head_.load(std::memory_order_relaxed);
			}
		}
	}

	size_t capacity() const { return n_; }

private:
	struct Slot {
		Slot() : seq(0), val() {}

		Slot(Slot&& other) noexcept
			: seq(0), val(std::move(other.val)) {
			size_t s = other.seq.load(std::memory_order_relaxed);
			seq.store(s, std::memory_order_relaxed);
		}
		Slot& operator=(Slot&& other) noexcept {
			if (this != &other) {
				size_t s = other.seq.load(std::memory_order_relaxed);
				seq.store(s, std::memory_order_relaxed);
				val = std::move(other.val);
			}
			return *this;
		}

		Slot(const Slot&) = delete;
		Slot& operator=(const Slot&) = delete;

		std::atomic<size_t> seq;
		T                   val;
	};

	std::vector<Slot>   slots_;
	size_t              n_{0}, mask_{0};
	std::atomic<size_t> head_{0}, tail_{0};
};


// ---------------------------------------
// Blocking wrapper: sleep only when needed
// ---------------------------------------
template<class T>
class BlockingQueue {
public:
	explicit BlockingQueue(size_t cap)
		: ring_(cap),
		  slots_(static_cast<int>(ring_.capacity())),
		  items_(0) {}

	void push(T v) {
		// acquire a free slot (fast path: try_acquire)
		if (!slots_.try_acquire()) {
			slots_.acquire();
		}

		// we logically own one slot now; enqueue must eventually succeed
		while (!ring_.try_enqueue(std::move(v))) {
			// contention in ring CAS, but capacity is guaranteed
			std::this_thread::yield();
		}

		// new item is now visible
		items_.release();
	}

	T pop() {
		// acquire an item
		if (!items_.try_acquire()) {
			items_.acquire();
		}

		T out;
		while (!ring_.try_dequeue(out)) {
			// someone else got there first; but we *know* an item exists
			std::this_thread::yield();
		}

		// freed one slot
		slots_.release();
		return out;
	}

	bool try_push(T v) {
		if (!slots_.try_acquire())
			return false;

		if (!ring_.try_enqueue(std::move(v))) {
			// logically took a slot but failed ring insert → undo
			slots_.release();
			return false;
		}
		items_.release();
		return true;
	}

	bool try_pop(T& out) {
		if (!items_.try_acquire())
			return false;

		if (!ring_.try_dequeue(out)) {
			// undo logical item take
			items_.release();
			return false;
		}
		slots_.release();
		return true;
	}

	size_t capacity() const { return ring_.capacity(); }

private:
	MpmcRing<T>          ring_;
	LightweightSemaphore slots_; // free capacity
	LightweightSemaphore items_; // items available
};

