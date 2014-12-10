/**
 * \file sym_iwt_2d.c
 *
 * 2D Biorthogonal Inverse Wavelet Transform
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

/*******************************************/
	void  sym_iwt_2d(sez_wcomp_par* S)
	{
        int    ii,jj,ik,jk,kk,ir,ic,jr,jc;
	int    qh,len,lev,ev,od,frst,scnd;
	float  temp, low_pass, hgh_pass;       
        float  qLS_1,qLS_2,qLS_3,qLS_4;
        float  qHS_0,qHS_1,qHS_2,qHS_3,qHS_4;
	float* aptr;

        qh  = (int) (S->nq)/2;

        qLS_1 = (S->qLS)[1];
        qLS_2 = (S->qLS)[2];
        qLS_3 = (S->qLS)[3];
        qLS_4 = (S->qLS)[4];
        qHS_0 = (S->qHS)[0];
        qHS_1 = (S->qHS)[1];
        qHS_2 = (S->qHS)[2];
        qHS_3 = (S->qHS)[3];
        qHS_4 = (S->qHS)[4];

        if (S->lev_r != 0) {
	/* multilevel row processing */  
        len  = (S->jr_0 - S->ir_0 + 1)/pow(2,S->lev_r);
        for (lev=1;lev<=S->lev_r;lev++){ 
          	for (ic=S->ic_0;ic <= S->jc_0;ic++) {
			frst = ic*(S->nr_ext)+(S->ir_0);
			scnd = frst+len;             
			/* interleave rows  */  
        		for (kk=0;kk<len;kk++){
				ik = frst+2*kk;
				jk = ik+1;
                		(S->Wou)[ik]=(S->Win)[frst+kk];
               			(S->Wou)[jk]=(S->Win)[scnd+kk];
          		}
			/* sym extension - rows */  
			scnd = frst+2*len-1;             
        		for (ir=1;ir <= qh;ir++) {
          			ii = frst-ir;
          			ik = frst+ir;
          			jj = scnd+ir;
          			jk = scnd-ir;
                      	        (S->Wou)[ii] = (S->Wou)[ik];
                        	(S->Wou)[jj] = (S->Wou)[jk];
			} 
			/* sym filtering - rows */  
                	for (ir=0;ir<2*len;ir=ir+2){
                        	jj = frst-qh+ir;
                                low_pass =((S->Wou)[jj+1]+(S->Wou)[jj+7])*qLS_1+((S->Wou)[jj+2]+(S->Wou)[jj+6])*qLS_2; 
                                low_pass+=((S->Wou)[jj+3]+(S->Wou)[jj+5])*qLS_3+(S->Wou)[jj+4]*qLS_4; 
                                hgh_pass =((S->Wou)[jj+1]+(S->Wou)[jj+9])*qHS_0+((S->Wou)[jj+2]+(S->Wou)[jj+8])*qHS_1; 
                                hgh_pass+=((S->Wou)[jj+3]+(S->Wou)[jj+7])*qHS_2+((S->Wou)[jj+4]+(S->Wou)[jj+6])*qHS_3; 
                       		hgh_pass+= (S->Wou)[jj+5]*qHS_4; 
				(S->Win)[frst+ir]   = low_pass;
				(S->Win)[frst+1+ir] = hgh_pass;
        	        }
                }
		len=2*len;
        }
        /* transposition */
          for (ir=S->ir_1;ir <= S->jr_1;ir++) {
          for (ic=S->ic_0;ic <= S->jc_0;ic++) {
           (S->Wou)[ir*(S->nc_ext)+ic] =  (S->Win)[ic*(S->nr_ext)+ir];
          }
          }
        aptr = S->Win;
        S->Win = S->Wou;
        S->Wou = aptr;
        }

        if (S->lev_c != 0) {
                                                                                                                             
	/* multilevel row processing */  
        len  = (S->jc_0 - S->ic_0 + 1)/pow(2,S->lev_c);
        for (lev=1;lev<=S->lev_c;lev++){ 
          	for (ir=S->ir_0;ir <= S->jr_0;ir++) {
			frst = ir*(S->nc_ext)+(S->ic_0);
			scnd = frst+len;             
			/* interleave rows  */  
        		for (kk=0;kk<len;kk++){
				ik = frst+2*kk;
				jk = ik+1;
                		(S->Wou)[ik]=(S->Win)[frst+kk];
               			(S->Wou)[jk]=(S->Win)[scnd+kk];
          		}
			/* sym extension - rows */  
			scnd = frst+2*len-1;             
        		for (ic=1;ic <= qh;ic++) {
          			ii = frst-ic;
          			ik = frst+ic;
          			jj = scnd+ic;
          			jk = scnd-ic;
                      	        (S->Wou)[ii] = (S->Wou)[ik];
                        	(S->Wou)[jj] = (S->Wou)[jk];
			} 
			/* sym filtering - rows */  
                	for (ic=0;ic<2*len;ic=ic+2){
				low_pass = .0; 
				hgh_pass = .0;
                        	jj = frst-qh+ic;
                                low_pass =((S->Wou)[jj+1]+(S->Wou)[jj+7])*qLS_1+((S->Wou)[jj+2]+(S->Wou)[jj+6])*qLS_2; 
                                low_pass+=((S->Wou)[jj+3]+(S->Wou)[jj+5])*qLS_3+(S->Wou)[jj+4]*qLS_4; 
                                hgh_pass =((S->Wou)[jj+1]+(S->Wou)[jj+9])*qHS_0+((S->Wou)[jj+2]+(S->Wou)[jj+8])*qHS_1; 
                                hgh_pass+=((S->Wou)[jj+3]+(S->Wou)[jj+7])*qHS_2+((S->Wou)[jj+4]+(S->Wou)[jj+6])*qHS_3; 
                       		hgh_pass+= (S->Wou)[jj+5]*qHS_4; 
				(S->Win)[frst+ic]   = low_pass;
				(S->Win)[frst+1+ic] = hgh_pass;
        	        }
                }
		len=2*len;
        }
        }
        }

