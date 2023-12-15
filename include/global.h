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

//global time
extern std::chrono::high_resolution_clock::time_point start;

// flags are global variables
extern short int tabular_flag;
extern int verbose_flag;
extern int maxThreads;

extern Log lg;
extern std::mutex mtx;
extern ThreadPool<std::function<bool()>> threadPool;
extern UserInput userInput;

#endif /* GLOBAL_H */
