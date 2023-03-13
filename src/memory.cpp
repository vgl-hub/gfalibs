
#include <sys/resource.h>

unsigned long int get_mem_usage(){
    
    struct rusage thisUsage;
    
    getrusage(RUSAGE_SELF, &thisUsage);
    return thisUsage.ru_maxrss;
    
}
