#include <stdint.h>
#include <cmath>
#include <atomic>

const char* memUnit[4] = {"B", "KB", "MB", "GB"};
std::atomic<int64_t> alloc(0), freed(0);

double get_mem_inuse(uint8_t unit){
    
    return (alloc - freed) / pow(1024, unit);
    
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
