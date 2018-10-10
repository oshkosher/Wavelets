#ifndef __NIX_TIMER__
#define __NIX_TIMER__

/*
   High-precision timer for Windows or Linux.

   Ed Karrels, ekarrels@scu.edu
*/

#ifdef _WIN32
#include <windows.h>
#endif
#include <errno.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <ctime>


// Just one method, NixTimer::time()
// Returns a high-precision relative time in seconds
class NixTimer {

#ifndef _WIN32

 public:
  static double time() {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t)) {
      fprintf(stderr, "Error calling clock_gettime: %s\n", strerror(errno));
      exit(1);
    }

    return t.tv_sec + 1e-9 * t.tv_nsec;
  }

#else // defined(_WIN32)

  LARGE_INTEGER frequency;
  static NixTimer nixTimerFrequencyObject;

 public:
  NixTimer() {
    if (!QueryPerformanceFrequency(&frequency))
      fprintf(stderr,"High-resolution performance counter not supported.",
	      true);
  }

  static double time() {
    LARGE_INTEGER counter;
    if (!QueryPerformanceCounter(&counter))
      fprintf(stderr, "error calling QueryPerformanceCounter");

    return counter.QuadPart / (double) nixTimerFrequencyObject.frequency.QuadPart;
  }
    
#endif
};    


#endif // __NIX_TIMER__
