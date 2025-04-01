#ifndef GLOBAL_H
#define GLOBAL_H

#include <mutex>
#include <chrono>
#include <queue>
#include <thread>
#include <functional>

#include "log.h"
#include "threadpool.h"
#include "bed.h"
#include "struct.h"
#include "bit-packing.h"

//global time
extern std::chrono::high_resolution_clock::time_point start;

//global variables
extern int verbose_flag;
extern Log lg;
extern std::vector<Log> logs; // log storage for verbose output. Each log in the vector comes from a separate job
extern int tabular_flag;

extern int maxThreads;
extern std::mutex mtx;
extern ThreadPool<std::function<bool()>> threadPool;

#endif /* GLOBAL_H */
