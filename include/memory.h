#ifndef MEMORY_H
#define MEMORY_H

extern const char* memUnit[4];
extern std::atomic<int64_t> alloc, freed;
extern uint64_t maxMem;
extern std::atomic<bool> freeMemory;

double get_mem_inuse(uint8_t unit);

double get_mem_usage(uint8_t unit);

double get_mem_total(uint8_t unit);

double convert_memory(int64_t value, uint8_t unit);

bool allocMemory(int64_t amount);

#endif //MEMORY
