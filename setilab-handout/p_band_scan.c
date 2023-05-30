#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <pthread.h>
#include <sched.h>
#include <assert.h>
#include <unistd.h> 
#include <math.h>
#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0
pthread_t* tid;
int num_processors;

typedef struct args {
    signal* sig;
    int band;
    double bandwidth;
    int filter_order;
    double* filter_coeffs;
    int threadCount;
    int num_bands;
    int bandsperthread;
    int leftovers;
    double* band_power;
  } args;


void* worker(void* arg){
  
  args* arguments = (args*)arg;
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(arguments->band % num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }
  int tidw = arguments->band;
  int bpt = arguments->bandsperthread;
  double bw = arguments->bandwidth;
  int fo = arguments->filter_order;
  int lo = arguments->leftovers;
  for(int band = tidw * bpt;
          band < (tidw) * bpt + lo;
          band += 1){
    int filter_index = band * (fo+1);
    generate_band_pass(arguments->sig->Fs,
                       band * bw + 0.0001, // keep within limits
                       (band + 1) * bw - 0.0001,
                       fo,
                       &(arguments->filter_coeffs)[filter_index]);
    hamming_window(arguments->filter_order,arguments->filter_coeffs);
    convolve_and_compute_power(arguments->sig->num_samples,
                               arguments->sig->data,
                               fo,
                               &(arguments->filter_coeffs)[filter_index],
                               &(arguments->band_power[band]));
   }
   pthread_exit(NULL); 
}


void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub,int num_threads) {

  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  double filter_coeffs[(filter_order + 1) * num_bands]; //Prevent race condition in filter gen
  double band_power[num_bands];

  tid = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
  int n_threads = num_threads;
  int bandsperthread;

  if (num_bands < num_threads){
    n_threads = num_bands;
    bandsperthread = 1;
  } else {
    bandsperthread = floor(num_bands/num_threads);
  }
  //printf("%d",n_threads);
 
 for (long i = 0; i < n_threads; i++) {
    args* threadArgs = malloc(sizeof(args));
    threadArgs->sig = sig;
    threadArgs->bandwidth = bandwidth;
    threadArgs->filter_order = filter_order;
    threadArgs->filter_coeffs = filter_coeffs;
    threadArgs->band_power = band_power;
    threadArgs->threadCount = n_threads;
    threadArgs->num_bands = num_bands;
    threadArgs->bandsperthread = bandsperthread;
    threadArgs->band = i;
    if (i == n_threads - 1){
      threadArgs->leftovers = bandsperthread + num_bands - bandsperthread*n_threads;
    } else {
      threadArgs->leftovers = bandsperthread;
    }
    int returncode = pthread_create(&(tid[i]),
                                    NULL,
                                    worker,
                                    (args*)threadArgs
                                    );
    if (returncode != 0) {
      perror("Failed to start thread");
      exit(-1);
    }
  }
  for (int i = 0; i < n_threads; i++) {
    int returncode = pthread_join(tid[i], NULL);
    if (returncode != 0) {
      perror("join failed");
      exit(-1);
    }
  }
  //free(threadArgs);
  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n"
         "User time        %lf seconds\n"
         "System time      %lf seconds\n"
         "Page faults      %ld\n"
         "Page swaps       %ld\n"
         "Blocks of I/O    %ld\n"
         "Signals caught   %ld\n"
         "Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  /*if (argc != 6) {
    usage();
    return -1;
  }*/

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  int num_threads  = atoi(argv[6]);
  num_processors     = atoi(argv[7]);
  if (num_processors > 32) {
    num_processors = 32;
  }
  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n"
         "file:     %s\n"
         "Fs:       %lf Hz\n"
         "order:    %d\n"
         "bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end,num_threads)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

