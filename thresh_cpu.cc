#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>    // std::sort
#include "thresh_cpu.h"

// Allocates absData on the heap and loads it with abs(data)
// Sorts absData.  Finds the index corresponding to the compRatio, and 
// reads that element to get the threshold. Sets min and max values that 
// remain above threshold  
float thresh_cpu(int len, float *data, float compRatio, float *maxVal)
{
	int loopCount = 10;
	int displayCount = 10;
	int count = len*len;	//Square matrix
	std::vector<float> absData;
	float maxV = -1;
	float minV = 1200;

	for (int idx = 0; idx < count; idx++ )
	{
		float fval = *(data+idx);
		float afval = (float)(abs(*(data+idx)));
		absData.push_back(afval);
		if ( loopCount >= 0)
		{
			std::cout << "idx: " << idx << " data[idx]:" << data[idx] << " absData[idx]: " << abs(data[idx]) << std::endl;
			loopCount--;   
		}
		if (afval > maxV) maxV = afval;
		if (afval < minV) minV = afval; 
	}
	
	int threshIdx = floor((compRatio) * count );
	std::cout << "threshIdx: " << threshIdx << " minV:" << minV << " maxV:" << maxV << std::endl;
	
	std::sort(absData.begin(),absData.end());
	std::cout << "myvector contains:";
	//for (std::vector<float>::iterator it=absData.begin(); it!=absData.end(); ++it)
	//	std::cout << ' ' << *it << std::endl;

	std::vector<float>::iterator it = absData.begin();
	it += (threshIdx-10);
	for (int idx = 0; idx < 20; idx++ )
	{
		std::cout << "idx:" << idx << " *it:" << *it << std::endl;
		it++;
	}

	it = absData.begin();
	it += threshIdx;	
	float threshold = *it;
	threshold = 0.90;
	
    std::cout << "Front:" << absData.front() << " Back:" << absData.back() << " threshold: " << threshold << std::endl;
	*maxVal = *(absData.rbegin()); // get the largest value in the data
        threshold += 1E-16;
	
	if (displayCount > 0 )
	{
		displayCount--;
		std::cout << "len:" << len << " maxVal:" << *maxVal << " threshIdx:" << threshIdx << " count:" << count << " threshold:" << threshold << std::endl;
	}

	return threshold;
}
