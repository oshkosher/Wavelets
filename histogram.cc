#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "data_io.h"

using namespace std;

struct Options {
  int columns, binCount;
  float trimFraction;  /* 0..1, not a percentage */
  const char *filename;
  int xOffset, yOffset, width, height;

  Options() {
    columns = 78;
    binCount = 24;
    trimFraction = 0.01f;
    filename = NULL;
    xOffset = 0;
    yOffset = 0;
    width = -1;
    height = -1;
  }
};


class FrequencyCounter {
  map<float,size_t> table;
  typedef map<float,size_t>::iterator TableIt;

public:
  void add(float f) {
    TableIt it = table.find(f);
    if (it == table.end()) {
      table[f] = 1;
    } else {
      it->second++;
    }
  }

  void getMode(float &value, size_t &count) {
    TableIt maxIt = table.begin();
    
    for (TableIt it = maxIt; it != table.end(); it++) {
      if (it->second > maxIt->second) maxIt = it;
    }

    value = maxIt->first;
    count = maxIt->second;
  }
};


void printHelp() {
  printf(
"\n"
"  histogram [options]  (read ASCII numbers from stdin)\n"
"  histogram [options] <input file> [<xOffset> <yOffset> <width> <height>]\n"
"    (auto-detect ASCII numbers or 2-d wavelet data file)\n"
"  options:\n"
"   -w <width>  Set the bar width, in characters\n"
"   -b <bins>   Set the number of bins\n"
"   -t <pct>    Set the percentage of outliers that are trimmed on each end\n"
"   -h          Output this help message\n"
"  If xOffset, yOffset, width, and height are supplied, then extract just\n"
"  a sub-rectangle of the given 2-d wavelet data file. If the input file\n"
"  is not a 2-d wavelet data file, output an error message.\n"
"\n");
  exit(1);
}


bool parseOptions(Options &opt, int argc, char **argv) {
  int argno = 1;

  while (argno < argc && argv[argno][0] == '-') {

    const char *arg = argv[argno];

    if (!strcmp(arg, "-w")) {
      argno++;
      if (argno == argc) printHelp();
      if (1 != sscanf(argv[argno], "%d", &opt.columns) ||
	  opt.columns <= 1) {
	printf("Invalid width: \"%s\"\n", argv[argno]);
	return false;
      }
    }
    
    else if (!strcmp(arg, "-b")) {
      argno++;
      if (argno == argc) printHelp();
      if (1 != sscanf(argv[argno], "%d", &opt.binCount) ||
	  opt.binCount < 1) {
	printf("Invalid count: \"%s\"\n", argv[argno]);
	return false;
      }
      
    }
    
    else if (!strcmp(arg, "-t")) {
      argno++;
      if (argno == argc) printHelp();
      if (1 != sscanf(argv[argno], "%f", &opt.trimFraction) ||
	  opt.trimFraction < 0 ||
	  opt.trimFraction >= 100) {
	printf("Invalid trim percentage: \"%s\"\n", argv[argno]);
	return false;
      }
      opt.trimFraction *= 0.01f;
    }

    else {
      printf("Invalid argument: \"%s\"\n", argv[argno]);
      printHelp();
    }

    argno++;
  }

  if (argno < argc) {
    // get the filename
    opt.filename = argv[argno++];

    // there are more arguments; must be xOffset..height
    if (argno < argc) {
      if (argc - argno != 4) printHelp();

      if (1 != sscanf(argv[argno], "%d", &opt.xOffset) ||
          1 != sscanf(argv[argno+1], "%d", &opt.yOffset) ||
          1 != sscanf(argv[argno+2], "%d", &opt.width) ||
          1 != sscanf(argv[argno+3], "%d", &opt.height)) {
        fprintf(stderr, "Cannot parse offsets\n");
        return false;
      }
    }
  }

  return true;
}


string makeStars(float count_f) {
  if (count_f == 0) {
    return "";
  } else {
    int count = (int) count_f;
    if (count == 0) {
      return "*";
    } else {
      return string(count, '*');
    }
  }
}


bool readFile(Options &opt, vector<float> &data) {

  data.clear();

  if (isValidDataFile(opt.filename)) {

    float *array;
    int width, height;

    if (!readDataFile(opt.filename, &array, &width, &height)) return false;
    
    if (opt.width < 0) opt.width = width - opt.xOffset;
    if (opt.height < 0) opt.height = height - opt.yOffset;

    if (opt.xOffset + opt.width > width ||
        opt.yOffset + opt.height > height) {
      fprintf(stderr, "Cannot extract %d..%d x %d..%d from %dx%d image.\n",
              opt.xOffset, opt.xOffset + opt.width - 1,
              opt.yOffset, opt.yOffset + opt.height - 1,
              width, height);
      delete[] array;
      return false;
    }

    data.reserve(opt.width * opt.height);

    // add each row of the data
    for (int y=opt.yOffset; y < opt.yOffset + opt.height; y++) {
      float *rowStart = array + y*width + opt.xOffset;
      data.insert(data.end(), rowStart, rowStart+opt.width);
    }

    delete[] array;
  }

  // plain text numeric input
  else {
    FILE *inf = stdin;

    if (opt.filename) {
      inf = fopen(opt.filename, "r");
      if (!inf) {
        fprintf(stderr, "Cannot read \"%s\"\n", opt.filename);
        return false;
      }
    }

    float f;
    while (true) {
      if (1 != fscanf(inf, "%f", &f)) break;
      data.push_back(f);
    }
    if (inf != stdin) fclose(inf);
  }

  return true;
}
  

int main(int argc, char **argv) {

  Options opt;
  vector<float> data;
  vector<int> binSizes;
  float min, max, trimmedMin, trimmedMax;
  size_t startIdx, endIdx;
  double sum = 0, sum2 = 0;

  if (!parseOptions(opt, argc, argv)) return 1;

  if (!readFile(opt, data)) return 1;

  if (data.size() == 0) {
    printf("No input.\n");
    return 0;
  }

  FrequencyCounter counter;

  // compute totals for average and standard deviation
  size_t n = data.size();
  for (size_t i=0; i < n; i++) {
    float f = data[i];
    sum += f;
    sum2 += f*f;
    counter.add(f);
  }

  sort(data.begin(), data.end());

  trimmedMin = min = data[0];
  trimmedMax = max = data[data.size()-1];

  if (min == max) {
    printf("Data is constant value: %g\n", min);
    return 0;
  }

  startIdx = 0;
  endIdx = data.size();

  // figure out how many entries to trim from either end
  if (opt.trimFraction > 0) {
    startIdx = (size_t) (opt.trimFraction * data.size());
    if (startIdx > 0) {
      endIdx -= startIdx;
      trimmedMin = data[startIdx];
      trimmedMax = data[data.size()-1-startIdx];
    }
  }

  float binWidth = (trimmedMax - trimmedMin) / opt.binCount;

  // scan the bins to compute all their sizes. The largest bin will then
  // be scaled to the screen size
  size_t binStartIdx = 0;
  float binStartValue = trimmedMin;
  float binEndValue = binStartValue + binWidth;
  size_t binEndIdx = lower_bound(data.begin(), data.end(), binEndValue)
    - data.begin();
				 
  size_t binSize = binEndIdx - binStartIdx;
  binSizes.push_back(binSize);
  binStartIdx = binEndIdx;
  binStartValue = binEndValue;
  size_t maxBinSize = binSize;

  for (int binId=1; binId < opt.binCount-1; binId++) {
    binEndValue = binStartValue + binWidth;
    binEndIdx = lower_bound(data.begin() + binStartIdx, data.end(),
			    binEndValue) - data.begin();
    binSize = binEndIdx - binStartIdx;
    binSizes.push_back(binSize);
    if (binSize > maxBinSize) maxBinSize = binSize;
    binStartIdx = binEndIdx;
    binStartValue = binEndValue;
  }
  binSize = data.size() - binStartIdx;
  binSizes.push_back(binSize);
  if (binSize > maxBinSize) maxBinSize = binSize;

  // scale to fit the largest bin
  float binSizeScale = (opt.columns - 13 + .5f) / maxBinSize;

  // output the bins
  binEndValue = trimmedMin + binWidth;
  string stars = makeStars(binSizes[0] * binSizeScale);
  printf(">= %-10.3g|%s\n", min, stars.c_str());
  binStartValue = binEndValue;

  for (int binId=1; binId < opt.binCount; binId++) {
    stars = makeStars(binSizes[binId] * binSizeScale);
    binEndValue = binStartValue + binWidth;
    printf(">= %-10.3g|%s\n", binStartValue, stars.c_str());
    binStartValue = binEndValue;
  }

  // compute the standard deviation
  double stddev = sqrt((sum2 - (sum*sum) / data.size()) / (data.size()-1));

  float mode;
  size_t modeCount;
  counter.getMode(mode, modeCount);

  // output some statistics
  printf("range %.7g .. %.7g\n", min, max);
  printf("count %llu, mode=%.7g (%llu)\n", (unsigned long long) data.size(),
	 mode, (unsigned long long) modeCount);
  printf("average %.7g\n", sum / data.size());
  printf("stddev %.7g\n", stddev);

  return 0;
}
