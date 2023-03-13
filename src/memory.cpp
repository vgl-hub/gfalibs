#ifndef _WIN32

#include <sys/resource.h>

#endif

unsigned long int get_mem_usage(){
    
#ifndef _WIN32
    
    struct rusage thisUsage;
    
    getrusage(RUSAGE_SELF, &thisUsage);
    return thisUsage.ru_maxrss;
    
#else
    
    return 0;
    
#endif
    
}
