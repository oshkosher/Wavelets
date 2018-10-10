/**
 * \file sym_fwt_2d.c
 *
 * 2D Biorthogonal Forward Wavelet Transform
 *
 * \author Sergio E. Zarantonello
 */
                                                                                
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "sez_wcomp.h"

#include <stdio.h>
void print_matrix(float *data, int rows, int cols, int rowWidth,
                  const char *name) {
  int row, col;
  if (name) puts(name);
  for (row=0; row < rows; row++) {
    for (col=0; col < cols; col++) {
      printf("%g, ", data[row*rowWidth + col]);
    }
    putchar('\n');
  }
}

/*******************************************/
	void  sym_fwt_2d(sez_wcomp_par* S)
	{
        int    ii,jj,ik,jk,kk,ir,ic,jr,jc;
	int    qh,len,lev,ev,od,frst,last;
	float  low_pass, hgh_pass;       
        float  qLA_0,qLA_1,qLA_2,qLA_3,qLA_4; 
        float  qHA_1,qHA_2,qHA_3,qHA_4; 
	float* aptr;

        qh  = (int) (S->nq)/2;

        qLA_0 = (S->qLA)[0];
        qLA_1 = (S->qLA)[1];
        qLA_2 = (S->qLA)[2];
        qLA_3 = (S->qLA)[3];
        qLA_4 = (S->qLA)[4];
        qHA_1 = (S->qHA)[1];
        qHA_2 = (S->qHA)[2];
        qHA_3 = (S->qHA)[3];
        qHA_4 = (S->qHA)[4];

	if (S->lev_c != 0){

          print_matrix(S->Win + S->ir_1*S->nc_ext + S->ic_0, 
                       S->jr_1 - S->ir_1 + 1,
                       S->jc_0 - S->ic_0 + 1,
                       S->nc_ext,
                       "before replication");

	/* left-right replication */
        for (ir=S->ir_1;ir <= S->jr_1;ir++) {
                kk = ir*(S->nc_ext)+(S->ic_1); 
        	for (ic=S->ic_0;ic < S->ic_1;ic++) {
          		ii = ir*(S->nc_ext)+ic;
                        (S->Win)[ii] = (S->Win)[kk];
		} 
                kk = ir*(S->nc_ext)+(S->jc_1); 
        	for (jc=S->jc_0;jc > S->jc_1;jc--) {
          		jj = ir*(S->nc_ext)+jc;
                        (S->Win)[jj] = (S->Win)[kk];
		} 
	}
        

        print_matrix(S->Win + S->ir_1*S->nc_ext + S->ic_0, 
                     S->jr_1 - S->ir_1 + 1,
                     S->jc_0 - S->ic_0 + 1,
                     S->nc_ext,
                     "before replication");

	/* multilevel row processing */  
        len  = S->jc_0 - S->ic_0 + 1;
        frst = S->ic_0;
        for (lev=1;lev<=S->lev_c;lev++){ 
                last = frst+len-1;
            	len  = (int) len/2;

		/* all rows */
        	for (ir=S->ir_1;ir <= S->jr_1;ir++) {

                  print_matrix(S->Win + ir*S->nc_ext + S->ic_0 - qh,
                               1, S->jc_0 - S->ic_0 + qh*2, S->nc_ext,
                               "before row processing");
                  

			/* symmetric extension */
        		for (ic=1;ic <= qh;ic++) {
          			ii = ir*(S->nc_ext)+frst-ic;
          			ik = ir*(S->nc_ext)+frst+ic;
          			jj = ir*(S->nc_ext)+last+ic;
          			jk = ir*(S->nc_ext)+last-ic;
                      	        (S->Win)[ii] = (S->Win)[ik];
                        	(S->Win)[jj] = (S->Win)[jk];
			} 

                        
                        print_matrix(S->Win + ir*S->nc_ext + S->ic_0 - qh,
                                     1, S->jc_0 - S->ic_0 + qh*2, S->nc_ext,
                                     "after symmetric extension");


			/* convolution with mirror filters */  
			ev = ir*(S->nc_ext)+qh;
			od = ev+len;
			kk = 0;
                	for (ic=0;ic<len;ic++){
                      		jj = ir*(S->nc_ext)+kk;
                       	low_pass = ( (S->Win)[jj]+(S->Win)[jj+8] )*qLA_0+( (S->Win)[jj+1]+(S->Win)[jj+7] )*qLA_1;
                       	low_pass +=( (S->Win)[jj+2]+(S->Win)[jj+6] )*qLA_2+( (S->Win)[jj+3]+(S->Win)[jj+5] )*qLA_3+(S->Win)[jj+4]*qLA_4; 
			hgh_pass = ((S->Win)[jj+2]+(S->Win)[jj+8])*qHA_1+((S->Win)[jj+3]+(S->Win)[jj+7])*qHA_2; 
                        hgh_pass +=((S->Win)[jj+4]+(S->Win)[jj+6])*qHA_3+(S->Win)[jj+5]*qHA_4;

                        printf("out[%d] = %d .. %d\n", ev, kk, kk+8);
                        printf("out[%d] = %d .. %d\n\n", od, kk+2, kk+8);
			(S->Wou)[ev] = low_pass;
			(S->Wou)[od] = hgh_pass;
			kk=kk+2;ev++;od++;
        	        }
			ev = ir*(S->nc_ext)+qh;
                	// for (ic=0;ic<len;ic++){
                	for (ic=0;ic<len;ic++){
				(S->Win)[ev]=(S->Wou)[ev];
				ev++;
			}
                        
                        print_matrix(S->Wou + ir*S->nc_ext + S->ic_0,
                                     1, S->jc_0 - S->ic_0 + 1, S->nc_ext,
                                     "Wou after transform");
                }
        }
	aptr = S->Win;
	S->Win = S->Wou;
	S->Wou = aptr;
	}

          print_matrix(S->Win + S->ir_1*S->nc_ext + S->ic_0, 
                       S->jr_1 - S->ir_1 + 1,
                       S->jc_0 - S->ic_0 + 1,
                       S->nc_ext,
                       "after row transform");

	if (S->lev_r != 0){

        /* transposition */ 
	  //for (ir=0;ir < S->nr_ext;ir++) {
          //for (ic=0;ic < S->nc_ext;ic++) {
	  for (ir=S->ir_1;ir <= S->jr_1;ir++) {
          for (ic=S->ic_0;ic <= S->jc_0;ic++) {
	   (S->Wou)[ic*(S->nr_ext)+ir] =  (S->Win)[ir*(S->nc_ext)+ic]; 
          } 
          } 
	/* left-right replication */
        for (ic=S->ic_0;ic <= S->jc_0;ic++){
                kk = ic*(S->nr_ext)+(S->ir_1); 
        	for (ir=S->ir_0;ir < S->ir_1;ir++) {
          		ii = ic*(S->nr_ext)+ir;
                        (S->Wou)[ii] = (S->Wou)[kk];
		} 
                kk = ic*(S->nr_ext)+(S->jr_1); 
        	for (jr=S->jr_0;jr > S->jr_1;jr--) {
          		jj = ic*(S->nr_ext)+jr;
                        (S->Wou)[jj] = (S->Wou)[kk];
		} 
	}
        /* multilevel row processing */
        len  = S->jr_0 - S->ir_0 + 1;
        frst = S->ir_0;
        for (lev=1;lev<=S->lev_r;lev++){
                last = frst+len-1;
                len  = (int) len/2;
                /* all rows */
                for (ic=S->ic_0;ic <= S->jc_0;ic++) {
                        /* symmetric extension */
                        for (ir=1;ir <= qh;ir++) {
                                ii = ic*(S->nr_ext)+frst-ir;
                                ik = ic*(S->nr_ext)+frst+ir;
                                jj = ic*(S->nr_ext)+last+ir;
                                jk = ic*(S->nr_ext)+last-ir;
                                (S->Wou)[ii] = (S->Wou)[ik];
                                (S->Wou)[jj] = (S->Wou)[jk];
                        }
                        /* convolution with mirror filters */
                        ev = ic*(S->nr_ext)+qh; 
                        od = ev+len;
                        kk = 0;
                        for (ir=0;ir<len;ir++){
                                jj = ic*(S->nr_ext)+kk;
                       	low_pass = ( (S->Wou)[jj]+(S->Wou)[jj+8] )*qLA_0+( (S->Wou)[jj+1]+(S->Wou)[jj+7] )*qLA_1;
                       	low_pass +=( (S->Wou)[jj+2]+(S->Wou)[jj+6] )*qLA_2+( (S->Wou)[jj+3]+(S->Wou)[jj+5] )*qLA_3+(S->Wou)[jj+4]*qLA_4; 
			hgh_pass = ((S->Wou)[jj+2]+(S->Wou)[jj+8])*qHA_1+((S->Wou)[jj+3]+(S->Wou)[jj+7])*qHA_2; 
                        hgh_pass +=((S->Wou)[jj+4]+(S->Wou)[jj+6])*qHA_3+(S->Wou)[jj+5]*qHA_4;

                        (S->Win)[ev] = low_pass;
                        (S->Win)[od] = hgh_pass;
                        kk=kk+2;ev++;od++;
                        }
			ev=ic*(S->nr_ext)+qh;
                        for (ir=0;ir<len;ir++){
                                (S->Wou)[ev]=(S->Win)[ev];
                                ev++;
                        }
                }
	}
	}
        }	
