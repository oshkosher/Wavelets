#ifndef __CUDA_TIMER_H__
#define __CUDA_TIMER_H__

/*
  CudaTimer - make it easier to use CUDA events to track the time
  spent during some operations.  Rather than creating a start event,
  recording it in the current stream, then creating an end event and
  recording it, then deallocating both, just make a CudaTime object
  and call start() and end().  If multiple calls to start() and end()
  are made, the time between each pair will be added up.
*/


#include <vector>
#include <string>
#include <cstdio>
#include "cuda.h"

using namespace std;

class CudaTimer {

  string name;
  vector<cudaEvent_t> events;
  bool isStarted;
  cudaStream_t currentStream;

  // add an event
  void add(cudaStream_t stream = 0) {
    cudaEvent_t event;
    CUCHECK(cudaEventCreate(&event));
    CUCHECK(cudaEventRecord(event, stream));
    events.push_back(event);
  }

 public:
  CudaTimer(const char *name_ = "") {
    name = name_;
    isStarted = false;
  }

  ~CudaTimer() {
    for (size_t i=0; i < events.size(); i++)
      cudaEventDestroy(events[i]);
  }

  const char *getName() {return name.c_str();}

  // start timer
  void start(cudaStream_t stream = 0) {
    if (isStarted) {
      fprintf(stderr, "ERROR: CudaTimer::start() called with timer already "
              "started.\n");
      return;
    }

    isStarted = true;
    currentStream = stream;
    add(stream);
  }

  // end timer
  void end(cudaStream_t stream = 0) {
    if (!isStarted) {
      fprintf(stderr, "ERROR: CudaTimer::end() called with timer not "
              "started.\n");
      return;
    }

    if (stream != currentStream) {
      fprintf(stderr, "ERROR: CudaTimer::end() called on a different stream "
              "than the one on which the timer was started.\n");
      return;
    }

    currentStream = 0;
    isStarted = false;
    add(stream);
  }

  // return the number of start/end pairs
  int count() {return events.size() / 2;}

  // Compute the total time in milliseconds.
  float time() {
    float total = 0;
    for (size_t i = 0; i < events.size(); i += 2) {
      float ms;
      CUCHECK(cudaEventElapsedTime(&ms, events[i], events[i+1]));
      total += ms;
    }
    return total;
  }

  // return the number of CUDA events
  int countEvents() {return events.size() / 2;}
  cudaEvent_t getEvent(size_t i) {return events[i];}
  cudaEvent_t getLastEvent() {return events[events.size()-1];}

  // return the time between event 2*i and 2*i+1
  float getEventTime(int i) {
    float ms;
    CUCHECK(cudaEventElapsedTime(&ms, events[i*2], events[i*2+1]));
    return ms;
  }
    
};

#endif // __CUDA_TIMER_H__
