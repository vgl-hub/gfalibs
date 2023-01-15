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
    short unsigned int running();
    unsigned int queueSize();
    void join();
    short unsigned int totalThreads();

friend class InSequences;
    
};

template<class T>
void ThreadPool<T>::threadLoop(int threadN) {
    
    T job;
    
    while (true) {
        
        {
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
            
            threadStates[threadN] = false;
            job = jobs.front();
            jobs.pop();
        }
        job();
        threadStates[threadN] = true;
#ifdef DEBUG
        std::cout<<"Thread "<<std::to_string(threadN)<<" done (thread state: "<<threadStates[threadN]<<")"<<std::endl;
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
    this->maxThreads = maxThreads-1;
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
short unsigned int ThreadPool<T>::running() {
    
    short unsigned int jobN = 0;
    
    for(bool isDone : threadStates) {
        if (!isDone)
            ++jobN;
    }
    
    return jobN;

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
short unsigned int ThreadPool<T>::totalThreads() {
    return maxThreads;
}

template<class T>
void jobWait(ThreadPool<T>& threadPool) {
    while (true) {
        
        if (threadPool.empty() && threadPool.jobsDone()) {
            
            lg.newlines(2);
            break;
            
        }
        lg.verbose("Jobs waiting/running: " + std::to_string(threadPool.queueSize()) + "/" + std::to_string(threadPool.running()), true);
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        
    }
}

#endif //THREADPOOL
