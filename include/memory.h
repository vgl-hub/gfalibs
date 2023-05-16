#ifndef MEMORY_H
#define MEMORY_H

extern const char* memUnit[4];
extern uint64_t freed;

double get_mem_usage(uint8_t unit);

double get_mem_total(uint8_t unit);

#endif //MEMORY
