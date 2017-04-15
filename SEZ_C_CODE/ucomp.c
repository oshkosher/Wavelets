/**
 * \file ucomp.c
 * 
 * 2D lossy+lossless compress/expand
 *
 *
 * \author Sergio Zarantonello
 */
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include "zlib_wrapper.h"
#include "sez_wcomp.h"
#include "sez_alloc.h"

void ucomp(float* dat_in, float* dat_ou,int* par_int)
{
    sez_wcomp_par* S;
    float v,f,max,min,*qwc,*qcc;
    int nr,nc,ns,nt,nrx0,ncx0,nrx1,ncx1;
    int i,j,ir,ic,is,frst,sv,jj;
    int wrds,shft,nn,na,nb,sbs_c,sbs_r,ind_min_wc;
    int len_wc,n,int_in,qlen,bytes_in,bytes_ou; 
    
    //nr=nx;nc=nz;ns=ny;
    nc  = par_int[0];
    nr  = par_int[1];
    ns  = par_int[2];
    qlen= par_int[3];
    wrds= par_int[5];

    S=sez_wcomp_init(nr,nc,1,1);

    ncx0=(S->jc_0)-(S->ic_0)+1;
    nrx0=(S->jr_0)-(S->ir_0)+1;
    ncx1=(S->jc_1)-(S->ic_1)+1;
    nrx1=(S->jr_1)-(S->ir_1)+1;

    len_wc=ncx0*nrx0;
    qwc = sez_float_vector(0,len_wc-1); 
    qcc = sez_float_vector(0,len_wc-1); 

    shft=0;jj=0;
    for (is=0;is<ns;is++){
    n=0;

      min      =       dat_in[shft]; 
      max      =       dat_in[shft+1]; 
      bytes_in = (int) dat_in[shft+2]; 

      //printf(" UCOMP %f\t %f\t %d\n",min,max,bytes_in);

      int_in   = (int) ceil((int) bytes_in/sizeof(int)  + 1);

      for (i=0;i<int_in;i++){qcc[i]= dat_in[i+shft+3];}

      shft = shft+int_in+3;
      bytes_ou = sizeof(int)*len_wc;

      zlib_uncompress(qwc,bytes_ou,qcc,bytes_in);

      if (max > min)
      {
        double f = (double) (max-min)/qlen;
        for (ic=S->ic_0;ic<=S->jc_0;ic++) {
          frst = ic*(S->nr_ext);
          for (ir=S->ir_0;ir<=S->jr_0;ir++) {
            v = (double) qwc[n++];
            (S->Win)[frst+ir] = v*f + min;
          }
        }
      }
      else
      {
        for (ic=S->ic_0;ic<=S->jc_0;ic++) {
          frst = ic*(S->nr_ext);
          for (ir=S->ir_0;ir<=S->jr_0;ir++) {
            (S->Win)[frst+ir] = min;
          }
        }
      }

      // inverse biorth
      sym_iwt_2d(S);

      for (i=S->ir_1;i<=S->jr_1;i++)
      {
        for (j=S->ic_1;j<=S->jc_1;j++)
        {
        dat_ou[jj++] = (S->Win)[(S->nc_ext)*i+j];
        }
      }
    }
    sez_free_vector(S->Win,0,(S->nc_ext)*(S->nr_ext)-1);
    sez_free_vector(S->Wou,0,(S->nc_ext)*(S->nr_ext)-1);
    sez_free_vector(S->qLA,0,S->nq-1);
    sez_free_vector(S->qLS,0,S->nq-1);
    sez_free_vector(S->qHA,0,S->nq-1);
    sez_free_vector(S->qHS,0,S->nq-1);
    free(S);
}
