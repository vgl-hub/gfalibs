#include <stdint.h>
#include <cmath>
#include <atomic>
#include <iostream>

#include "global.h"

const char* memUnit[4] = {"B", "KB", "MB", "GB"};
std::atomic<int64_t> alloc(0), freed(0);
std::atomic<bool> freeMemory(false);
uint64_t maxMem = 0;

double get_mem_inuse(uint8_t unit){
    
    int64_t inUse = (alloc - freed);
    
    if (inUse > maxMem * pow(1024, unit)) {
        {
            std::lock_guard<std::mutex> lck(mtx);
            freeMemory = true;
        }
        threadPool.notify_all();
    }else if (inUse < maxMem * 0.1 && freeMemory){
        {
            std::lock_guard<std::mutex> lck(mtx);
            freeMemory = false;
        }
        threadPool.notify_all();
    }
    
    return (alloc - freed) / pow(1024, unit);
    
}

double convert_memory(int64_t value, uint8_t unit) {
    
    return value / pow(1024, unit);
    
}

#ifdef _WIN32

#include <windows.h>

double get_mem_usage(uint8_t unit){
    return 0;
}

double get_mem_total(uint8_t unit){
    
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
    
}

#else

#include <sys/resource.h>
#include <unistd.h>

double get_mem_usage(uint8_t unit){
    
    struct rusage thisUsage;
    getrusage(RUSAGE_SELF, &thisUsage);
    
#ifdef __linux__
    
    return thisUsage.ru_maxrss / pow(1024, unit - 1);

#else
    
    return thisUsage.ru_maxrss / pow(1024, unit);

#endif
    
}

double get_mem_total(uint8_t unit){
    
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    
    return pages * page_size / pow(1024, unit);
    
}

#endif

bool allocMemory(int64_t amount) {
    
    while (get_mem_inuse(3) + convert_memory(amount, 3) > maxMem){}
    alloc += amount;
    
    return true;
    
}
