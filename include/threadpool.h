#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <iostream>
#include <condition_variable>
#include <future>
#include "memory.h"

#include "parallel-hashmap/phmap.h"

extern Log lg;
extern std::vector<Log> logs;
extern std::mutex mtx;

template<class T>
struct JobWrapper {
    uint32_t jid;
    T job;
    
};

template<class T>
class ThreadPool {
private:
    int maxThreads;
    std::vector<std::thread> threads;
    std::vector<bool> threadStates;
    std::queue<JobWrapper<T>> jobs;
    std::mutex queueMutex;
    std::condition_variable workersMtxCondition, masterMtxCondition;
    bool done = false;
    uint32_t uid = 1;
    std::chrono::high_resolution_clock::time_point past;

    void threadLoop(int threadN);

public:
    phmap::flat_hash_map<uint32_t, bool> queueJids;
    
    void init(int maxThreads);
    uint32_t queueJob(const T& job);
    std::vector<uint32_t> queueJobs(const std::vector<T> &jobs);
    bool empty();
    bool jobsDone();
    short unsigned int running();
    unsigned int queueSize();
    void join();
    uint32_t totalThreads();
    void execJob();
    void status();
	std::mutex& getMutex();
    std::condition_variable& getWorkersMtxCondition();
	std::condition_variable& getMasterMtxCondition();
    void notify_all();
    
};

template<class T>
void ThreadPool<T>::threadLoop(int threadN) {
    
    T job;
    uint32_t jid = 0;
    
    while (true) {
        
        {
            std::unique_lock<std::mutex> lock(queueMutex);
#ifdef DEBUG
            std::cout<<"Thread "<<std::to_string(threadN)<<" waiting"<<std::endl;
#endif
            threadStates[threadN] = true;
            
			workersMtxCondition.wait(lock, [this] {
                return !jobs.empty() || done;
            });
			masterMtxCondition.notify_one();
            if (done)
                return;
        }
        while (true) {
			{
				std::lock_guard<std::mutex> lock(queueMutex);
				queueJids[jid] = false; // set the job as executed
				if (jobs.empty()) { // return to wait if no more jobs available
					masterMtxCondition.notify_one();
					break;
				}
                threadStates[threadN] = false; // thread unavailable
                JobWrapper<T> jobWrapper = jobs.front(); // get the job
                job = jobWrapper.job;
                jid = jobWrapper.jid;
                jobs.pop();
            }
            job(); // execute the job
#ifdef DEBUG
            std::cout<<"Thread "<<std::to_string(threadN)<<" done (thread state: "<<threadStates[threadN]<<")"<<std::endl;
#endif
			masterMtxCondition.notify_one();
        }
    }
}

template<class T>
void ThreadPool<T>::init(int maxThreads) {
    
    if(maxThreads == 0) maxThreads = std::thread::hardware_concurrency();
    if(maxThreads == 0) maxThreads = 1;
    threads.resize(maxThreads);
    threadStates.resize(maxThreads);
    
    lg.verbose("Generating threadpool with " + std::to_string(maxThreads) + " threads");
    
    for(int i=0; i<maxThreads; ++i) {
        threads[i] = std::thread(&ThreadPool::threadLoop, this, i);
        threadStates[i] = true;
    }
    this->maxThreads = maxThreads;
    done = false;
}

template<class T>
uint32_t ThreadPool<T>::queueJob(const T& job) {
    
    uint32_t jid;
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        
        jid = ++uid;
        JobWrapper<T> jobWrapper{jid, job};
        jobs.push(jobWrapper);
        queueJids[jid] = true;
    }
	workersMtxCondition.notify_one();
    
    return jid;
}

template<class T>
std::vector<uint32_t> ThreadPool<T>::queueJobs(const std::vector<T> &newJobs) {
    
    std::vector<uint32_t> jids;
    uint32_t jid;
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        
        for (const T &job : newJobs) {
            
            jid = ++uid;
            JobWrapper<T> jobWrapper{jid, job};
            jobs.push(jobWrapper);
            queueJids[jid] = true;
            jids.push_back(jid);
            
        }
    }
	workersMtxCondition.notify_all();
    
    return jids;
}

template<class T>
bool ThreadPool<T>::empty() {
    
    if (!jobs.empty())
        workersMtxCondition.notify_all();
    
    return jobs.empty();
    
}

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
        std::lock_guard<std::mutex> lock(queueMutex);
        done = true;
    }
    workersMtxCondition.notify_all();
    for(std::thread& activeThread : threads) {
        activeThread.join();
    }
    threads.clear();
}

template<class T>
uint32_t ThreadPool<T>::totalThreads() {
    return maxThreads;
}

template<class T>
void ThreadPool<T>::execJob() {
    
    T job;
    uint32_t jid = 0;
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        
        queueJids[jid] = false;
        
        if (empty())
            return;
            
        JobWrapper<T> jobWrapper = jobs.front();
        job = jobWrapper.job;
        jid = jobWrapper.jid;
        jobs.pop();

    }
    job();

}

template<class T>
void ThreadPool<T>::status() {
    
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
    
    if (elapsed.count() > 0.1) {
        lg.verbose("Jobs waiting/running: " + std::to_string(queueSize()) + "/" + std::to_string(running()) + " memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
        
        past = std::chrono::high_resolution_clock::now();
    }
    
}

template<class T>
std::mutex& ThreadPool<T>::getMutex() {
	return queueMutex;
}

template<class T>
std::condition_variable& ThreadPool<T>::getWorkersMtxCondition() {
    return workersMtxCondition;
}

template<class T>
std::condition_variable& ThreadPool<T>::getMasterMtxCondition() {
	return masterMtxCondition;
}

template<class T>
void ThreadPool<T>::notify_all() {
    workersMtxCondition.notify_all();
}

inline void flushLogs() {
    
    std::unique_lock<std::mutex> lck(mtx);
    for (auto it = logs.begin(); it != logs.end(); it++) {
        it->print();
        logs.erase(it--);
    }
}

template<class T>
void jobWait(ThreadPool<T>& threadPool) {
    
    uint32_t jobNumber = threadPool.queueSize();
    
    while (true) {
        
        if (jobNumber > threadPool.queueSize()) {
            flushLogs();
            jobNumber = threadPool.queueSize();
        }
        threadPool.status();
        
		if (!threadPool.empty() || !threadPool.jobsDone()) {
			std::mutex &queueMutex = threadPool.getMutex();
			std::condition_variable& masterMtxCondition = threadPool.getMasterMtxCondition();
			std::unique_lock<std::mutex> lock(queueMutex);
			masterMtxCondition.wait(lock, [&threadPool] {
				return !threadPool.queueSize();
			});
		}else if (threadPool.empty() && threadPool.jobsDone()) {
			flushLogs();
			lg.verbose("\n", true);
			break;
		}
    }
}

template<class T>
void jobWait(ThreadPool<T>& threadPool, std::vector<uint32_t>& dependencies) {
    
    bool end = false;
    phmap::flat_hash_map<uint32_t, bool>::const_iterator got;
    uint32_t jobNumber = threadPool.queueSize();
    
    while (true) {

        threadPool.status();
        
        if (jobNumber > threadPool.queueSize()) {
            flushLogs();
            jobNumber = threadPool.queueSize();
        }
        for (uint32_t dependency : dependencies) {
            
            got = threadPool.queueJids.find(dependency);
            
            if (got == threadPool.queueJids.end()) {
                std::cout<<"Error: job dependency not found (id: "<<dependency<<")."<<std::endl;
                exit(0);
            }
            
            if (got->second) {
                end = false;
                break; // check this doesn't look right, i.e. case where dep is found
            }else{
                end = true;
            }
            
        }
        
        if (end == true) {
            lg.verbose("\n", true);
            dependencies.clear();
            break;
        }
        
		if (!threadPool.empty() || !threadPool.jobsDone()) {
			std::mutex &queueMutex = threadPool.getMutex();
			std::condition_variable& masterMtxCondition = threadPool.getMasterMtxCondition();
			std::unique_lock<std::mutex> lock(queueMutex);
			masterMtxCondition.wait(lock, [&threadPool] {
				return !threadPool.queueSize();
			});
		}else if (threadPool.empty() && threadPool.jobsDone()) {
			flushLogs();
			lg.verbose("\n", true);
			break;
		}
    }
}

#endif //THREADPOOL
