#ifndef THREADPOOL
#define THREADPOOL

#include <condition_variable>

extern Log lg;

template<class T>
class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::vector<bool> threadStates;
    std::queue<T> jobs;
    std::mutex queueMutex;
    std::condition_variable mutexCondition;
    bool done = false;

    void threadLoop(int threadN);

public:
    void init(int maxThreads);
    void queueJob(const T& job);
    bool empty();
    bool jobsDone();
    unsigned int queueSize();
    void join();

friend class InSequences;
    
};

template<class T>
void ThreadPool<T>::threadLoop(int threadN) {
    
    while (true) {
        T job;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
#ifdef DEBUG
            lg.verbose("Thread " + std::to_string(threadN) + " waiting");
#endif
            
            mutexCondition.wait(lock, [this] {
                return !jobs.empty() || done;
            });
            if (done) {
                return;
            }
            job = jobs.front();
            jobs.pop();
            
        }
        threadStates[threadN] = false;
        threadStates[threadN] = job();
#ifdef DEBUG
        lg.verbose("Thread " + std::to_string(threadN) + " done");
#endif

    }
}

template<class T>
void ThreadPool<T>::init(int maxThreads) {
    
    if(maxThreads == 0) maxThreads = std::thread::hardware_concurrency();
    if(maxThreads == 0 || maxThreads == 1) maxThreads = 2;
    threads.resize(maxThreads-1);
    threadStates.resize(maxThreads-1);
    
    lg.verbose("Generating threadpool with " + std::to_string(maxThreads-1) + " threads");
    
    for(int i=0; i<maxThreads-1; ++i) {
        threads[i] = std::thread(&ThreadPool::threadLoop, this, i);
        threadStates[i] = true;
    }
    this->maxThreads = maxThreads;
    done = false;
}

template<class T>
void ThreadPool<T>::queueJob(const T& job) {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        jobs.push(job);
    }
    mutexCondition.notify_one();
}

template<class T>
bool ThreadPool<T>::empty() {return jobs.empty();}

template<class T>
unsigned int ThreadPool<T>::queueSize() {return jobs.size();}

template<class T>
bool ThreadPool<T>::jobsDone() {
    
    for(bool done : threadStates) {
        if (!done)
//#ifdef DEBUG
            lg.verbose(done == 1 ? "done" : "not done");
//#endif
            return false;
    }
    
    return true;

}

template<class T>
void ThreadPool<T>::join() {
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        done = true;
    }
    mutexCondition.notify_all();
    for(std::thread& activeThread : threads) {
        activeThread.join();
    }
    threads.clear();
}

template<class T>
void jobWait(ThreadPool<T>& threadPool) {
    while (true) {
        
        if (threadPool.empty() && threadPool.jobsDone()) {break;}
        lg.verbose("Remaining jobs: " + std::to_string(threadPool.queueSize()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }
}

#endif //THREADPOOL
