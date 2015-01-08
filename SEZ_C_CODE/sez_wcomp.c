/**
 * \file sez_wcomp.c
 *
 * 2D Biorthogonal Wavelet Transform 
 *
 *
 * \author Sergio E. Zarantonello
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
                                                                                
#include "sez_wcomp.h"
#include "sez_alloc.h"

sez_wcomp_par* sez_wcomp_init(int nr, int nc, int lev_r, int lev_c)
{
sez_wcomp_par* S;
int nq_half;

float*  qA;
float*  qS;

/* Allocate structure */ 
   S = malloc(sizeof(sez_wcomp_par));

/* Input parameters */
    S->nc =nc;
    S->nr =nr;

/* Ajustable parameters */
    S->lev_r  = lev_r; /* processing steps for rows */
    S->lev_c  = lev_c; /* processing steps for columns */

/* Hard-wired filter length */
    S->nq = 9;

    /* cache the 'radius' of the filter, half its length */
    nq_half = (int) (S->nq)/2;

/*  Processing levels & offsets */

    // set up column indices
    if (nc >= S->nq && S->lev_c != 0){ 
        // if given -1, do the maximum number of levels
        if (S->lev_c==-1){
    	  S->lev_c = (int) floor( log( ((float) nc))/log(2) - 2);
        }

        // number of columns after expansion
        S->nc_ext = ceilf(nc/pow(2,S->lev_c))*pow(2,S->lev_c) + 2*nq_half; 

        // set up mirrored data
        S->ic_0 = nq_half;                        // expanded data start
        S->ic_1 = (int) ((S->nc_ext)-(S->nc))/2;  // initial data start
        S->jc_0 = (S->nc_ext) - nq_half - 1;      // expanded data end
        S->jc_1 = (S->nc) + (S->ic_1) -1;         // initial data end

    }
    else {
	S->lev_c=0;
        S->nc_ext = nc; 
        S->ic_0 = 0; 
        S->ic_1 = 0;
        S->jc_0 = nc-1; 
        S->jc_1 = nc-1; 
    }

    // set up row indices
    if (nr >= S->nq && S->lev_r != 0){ 
        if (S->lev_r==-1){
    	  S->lev_r =  (int) floor( log( ((float) nr))/log(2) - 2 );
          if (S->lev_r>4){S->lev_r=4;}
        }
        S->nr_ext = ceil(nr/pow(2,S->lev_r))*pow(2,S->lev_r) + 2*nq_half; 
        S->ir_0 = nq_half; 
        S->ir_1 = (int) ((S->nr_ext)-(S->nr))/2; 
        S->jr_0 = (S->nr_ext) - nq_half - 1; 
        S->jr_1 = (S->nr) + (S->ir_1) -1; 
    }
    else {
	S->lev_r=0;
        S->nr_ext = nr; 
        S->ir_0 = 0; 
        S->ir_1 = 0;
        S->jr_0 = nr-1; 
        S->jr_1 = nr-1; 
    }

/*  Allocate filters */
    S->qLA = sez_float_vector(0,8); 
    S->qHA = sez_float_vector(0,8); 
    S->qLS = sez_float_vector(0,8); 
    S->qHS = sez_float_vector(0,8); 

/* Hard-wired analysis and synthesis filters */
    qA = sez_float_vector(0,4); 
    qS = sez_float_vector(0,4); 

    qA[0] =  .037828455506995;
    qA[1] = -.02384946501938;
    qA[2] = -.110624404418420;
    qA[3] =  .377402855612650;
    qA[4] =  .85269867900940;

    qS[0] =  0;
    qS[1] = -.064538882628938;
    qS[2] = -.040689417609558;
    qS[3] =  .418092273222210;
    qS[4] =  .788485616405660;

/*  Filters for forward transform */
    (S->qLA)[0] = (S->qLA)[8] =  qA[0];
    (S->qLA)[1] = (S->qLA)[7] =  qA[1];
    (S->qLA)[2] = (S->qLA)[6] =  qA[2];
    (S->qLA)[3] = (S->qLA)[5] =  qA[3];
    (S->qLA)[4]               =  qA[4];

    (S->qHA)[0] = (S->qHA)[8] = -qS[0];
    (S->qHA)[1] = (S->qHA)[7] =  qS[1];
    (S->qHA)[2] = (S->qHA)[6] = -qS[2];
    (S->qHA)[3] = (S->qHA)[5] =  qS[3];
    (S->qHA)[4]               = -qS[4];

/*  Filters for inverse transform */
    (S->qLS)[0] = (S->qLS)[8] =  qS[0];
    (S->qLS)[1] = (S->qLS)[7] =  qA[1];
    (S->qLS)[2] = (S->qLS)[6] =  qS[2];
    (S->qLS)[3] = (S->qLS)[5] =  qA[3];
    (S->qLS)[4]               =  qS[4];

    (S->qHS)[0] = (S->qHS)[8] = -qA[0];
    (S->qHS)[1] = (S->qHS)[7] =  qS[1];
    (S->qHS)[2] = (S->qHS)[6] = -qA[2];
    (S->qHS)[3] = (S->qHS)[5] =  qS[3];
    (S->qHS)[4]               = -qA[4];

/* Allocate arrays for input/output */
   int length = S->nr_ext*S->nc_ext;
   S->Win = sez_float_vector(0,length-1);
   S->Wou = sez_float_vector(0,length-1);

   free(qA);free(qS);
   return S;
}
