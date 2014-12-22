 /******  sez_alloc.c
 * Basic memory allocation routines.
 ******/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
/****** error handler ******/
void nrerror(error_text)
        char error_text[];
        {
                void exit();
                fprintf(stderr,"run-time error\n");
                fprintf(stderr,"%s\n",error_text);
                fprintf(stderr,"exiting to system\n");
        }
/***** allocate matrix *****/
float **sez_float_matrix(int nrl, int nrh, int ncl, int nch)
        {
          int i;
          float **m;
          void nrerror();
          m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
          if (!m) nrerror("allocation failure 1 in matrix()");
          m -= nrl;
          for(i=nrl;i<=nrh;i++) {
                m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
                if(!m[i]) nrerror("allocation failure 2 in matrix()");
                m[i] -= ncl;
          }
          return m;
        }
/***** allocate vector *****/
float *sez_float_vector(int nl, int nh)
	{
	float *v;
	v =(float *) calloc(nh-nl+1, sizeof(float));
	/*v =(float *) malloc((unsigned) (nh-nl+1)*sizeof(float));*/
	if(!v) nrerror("allocation failure in sez_float_vector()");
	return v-nl;
	}
/***** allocate int vector *****/
int *sez_int_vector(int nl, int nh)
	{
	int *v;
	v =(int *) malloc((unsigned) (nh-nl+1)*sizeof(int));
	if(!v) nrerror("allocation failure in sez_int_vector()");
	return v-nl;
	}
/***** allocate short vector *****/
short *sez_short_vector(int nl, int nh)
	{
	short *v;
	v =(short *) malloc((unsigned) (nh-nl+1)*sizeof(short));
	if(!v) nrerror("allocation failure in sez_short_vector()");
	return v-nl;
	}
/****** free matrix ******/
void sez_free_matrix(float** m, int nrl, int nrh, int ncl, int nch)
	{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
	}
/****** free float vector ******/
void sez_free_vector(float* v, int nl, int nh)
	{
	free((char*) v+nl);
	}
/****** free int vector ******/
void sez_free_int(int* v, int nl, int nh)
	{
	free((char*) v+nl);
	}
/****** free short vector ******/
void sez_free_short(short* v, int nl, int nh)
	{
	free((char*) v+nl);
	}


