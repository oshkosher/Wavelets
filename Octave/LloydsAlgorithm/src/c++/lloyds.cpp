/*****************************************************************************/
/* C++ implementation of Lloyd's quantization algorithm                      */
/* lloyd(points, psize, codebook, csize, stop_criteria, partition,           */
/*		 dist, reldist, groups)                                              */
/*		                                                                     */
/*  - IN     - points:        values to be quantized                         */
/*  - IN     - psize:         size of the points array                       */
/*  - IN/OUT - codebook:      initial codebook and final codebook            */
/*  - IN     - csize:         size of the codebook array                     */
/*  - OUT    - partition:     array of csize-1 values with lloyd's partition */
/*  - OUT    - dist:          difference between each point and its          */
/*                            quantized value                                */
/*  - OUT    - reldist:       difference between consecutive dist. Used to   */
/*                            finish the iterating process                   */
/*  - OUT    - groups:        assign each point its codebook index           */
/*  - IN     - stop_criteria: typically 10e-7                                */
/*****************************************************************************/

void lloyd(float *points, unsigned int psize,float *codebook, unsigned int csize, float *partition, float &dist, float &reldist, unsigned int *groups, float stop_criteria = 10e-7) {
	float   *sum   = (float*)calloc(csize,sizeof(float));
	unsigned *count = (unsigned*)calloc(csize,sizeof(unsigned));

	// Calculate initial table
	for(int i = 0; i < csize-1; i++) {
		partition[i] = (codebook[i] + codebook[i+1]) / 2;
	}
	
	float pmax = points[0];
	float pmin = points[0];
	dist = 0;

	// Assign points to codebook (for codebook update)
	for(int i = 0; i < psize; i++) {
		float auxp = points[i];
		pmax = (auxp > pmax) ? auxp : pmax;
		pmin = (auxp < pmin) ? auxp : pmin;
		float d = abs(auxp - codebook[0]);
		float min = d;
		unsigned int idx = 0;
		for(int j = 1; j < csize; j++) {
			d = abs(auxp - codebook[j]);
			
			if (d < min) {
				idx = j;
				min = d;
			}
		}
		groups[i] = idx;
		sum[idx] += auxp;
		count[idx]++;
		dist += codebook[idx] - auxp;
	}
	dist /= psize;
	reldist = abs(dist);

	while(reldist > stop_criteria) {
	printf("reldist: %f(%f)\n",reldist,stop_criteria);

		// Update codebook
		for(int i = 0; i < csize; i++) {
			printf("sum: %f(%d)\n",sum[i],count[i]);
			if(count[i] != 0) {
				codebook[i] = sum[i] / count[i];
				count[i] = 0;
			} else if (i == 0) {
				codebook[i] = (partition[0] + pmin) / 2;
			} else if (i == csize-1) {
				codebook[i] = (partition[i-1] + pmax) / 2;
			} else {
				codebook[i] = (partition[i-1] + partition[i]) / 2;
			}
			sum[i] = 0;
		}

		// points in between codebooks (for table update)
		for(int i = 0; i < psize; i++) {
			for(int j = 0; j < csize-1; j++) {
				float auxp = points[i];
				if((auxp >= codebook[j])
				&& (auxp < codebook[j+1])) {
					count[j]++;
					sum[j] += auxp;
					break;
				}
			}
		}

		// Update table
		for(int i = 0; i < csize-1; i++) {
			if(count[i] != 0) {
				partition[i] = sum[i] / count[i];
				count[i] = 0;
			} else {
				partition[i] = (codebook[i] + codebook[i+1]) / 2;
			}
			sum[i] = 0;
		}

		reldist = dist;
		dist = 0;

		// Assign points to codebook (for codebook update)
		for(int i = 0; i < psize; i++) {
			float auxp = points[i];
			float d = abs(auxp - codebook[0]);
			float min = d;
			unsigned int idx = 0;
			for(int j = 1; j < csize; j++) {
				d = abs(auxp - codebook[j]);
				if (d < min) {
					idx = j;
					min = d;
				}
			}
			groups[i] = idx;
			sum[idx] += auxp;
			count[idx]++;
			dist += codebook[idx] - auxp;
		}
		dist /= psize;
		reldist = abs(reldist - dist);
	}  // END WHILE
}
