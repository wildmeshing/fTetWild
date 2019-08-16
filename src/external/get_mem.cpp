#include <floattetwild/get_mem.h>

extern "C" size_t getCurrentRSS();
extern "C" size_t getPeakRSS();

size_t floatTetWild::get_mem(){
    return getCurrentRSS()/(1024*1024);
}

size_t floatTetWild::get_peak_mem(){
    return getPeakRSS()/(1024*1024);
}