#ifndef THREADPOOL
#define THREADPOOL

#include <iostream>
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
    
};

template<class T>
void ThreadPool<T>::threadLoop(int threadN) {
    
    while (true) {
        
        threadStates[threadN] = true;

        std::unique_lock<std::mutex> lock(queueMutex);
#ifdef DEBUG
        std::cout<<"Thread "<<std::to_string(threadN)<<" waiting"<<std::endl;
#endif
        mutexCondition.wait(lock, [this] {
            return !jobs.empty() || done;
        });
        if (done) {
            return;
        }
        
        threadStates[threadN] = true;
        if(!jobs.empty()) {
            threadStates[threadN] = jobs.front()();
            jobs.pop();
        }
#ifdef DEBUG
        std::cout<<"Thread "<<std::to_string(threadN)<<" done"<<std::endl;
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
    
    for(bool isDone : threadStates) {
#ifdef DEBUG
        std::cout<<(isDone == true ? "done" : "not done")<<std::endl;
#endif
        if (!isDone)
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
