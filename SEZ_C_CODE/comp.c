/**
 * \file comp.c
 * 
 * 3D (by2D slice) lossy+lossless compress/expand
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

void comp(float* dat_in,int* par_int) 
/*
dat_in: input 1D array of dimensions nz*nx*ny.   
*/
{
   sez_wcomp_par* S;
   float v,min,max,cmp_factor,*qwc;
   int nr,nc,ns,nrx0,ncx0,nrx1,ncx1;
   int wrds,i,ir,ic,is,frst,sv;
   int n_in,n_ou,na,nb,sbs_c,sbs_r,ind_min_wc;
   int len_wc,qlen,bytes_in,int_ou; 
   size_t bytes_ou;
    
   nc  = par_int[0];
   nr  = par_int[1];
   ns  = par_int[2];
   qlen= par_int[3];

   // structure containing all the parameters
   // -1 means do the maximum number of wavelet steps
   S=sez_wcomp_init(nr,nc,-1,-1);

   // define an offset into the data to use different indices

   // number of columns and rows after expansion
   // i = first element, j = last element
   ncx0=(S->jc_0)-(S->ic_0)+1;
   nrx0=(S->jr_0)-(S->ir_0)+1;

   // actual number of elements in a column and row
   ncx1=nc;nrx1=nr;

   // length of the buffer in 1 dimension
   len_wc=ncx0*nrx0;

   // allocate temporary storage array
   qwc = sez_float_vector(0,len_wc-1); 

   wrds=0;n_in=0;n_ou=0;

   // ns == the Z dimension of the cube
   // each 2-d cross section is processed independently
   for (is=0;is<ns;is++)
   { 
      nb=0;  // index into the output array

      // copy from dat_in[] to S->Win[]
      for (ir=S->ir_1;ir<=S->jr_1;ir++)  // rows
      {
        for (ic=S->ic_1;ic<=S->jc_1;ic++)  // columns
        {
        (S->Win)[(S->nc_ext)*ir+ic]=dat_in[n_in++];
        }
      }

      sym_fwt_2d(S);
                                                                                
      // find min and max wavelet coeffs
      min = 1.0E20; max = -1.0E20;
      for (ic=S->ic_0;ic<=S->jc_0;ic++) {
        frst=ic*(S->nr_ext);
        for (ir=S->ir_0;ir<=S->jr_0;ir++) {
          v = (S->Win)[frst+ir];
          if(v>=max) max=v; if(v<=min) min=v;
        }
      }

      // quantize all wavelet coeffs (uniform quantization)
      // qlen is the range of quantized values
      if (min< max){
	 double f = (double) qlen / (max-min);
         for (ic=S->ic_0;ic<=S->jc_0;ic++) 
         {
           frst=ic*(S->nr_ext);
           for (ir=S->ir_0;ir<=S->jr_0;ir++) 
           {
             float v=(S->Win)[frst+ir];
             qwc[nb++] =  ceil((v-min)*f-.5);
           }
         }
      } 
      else 
      {
         // if min==max, just fill with zeros
         // XXX little bug--this should loop over all the data
         qwc[nb++] = 0.;
      }
      
      bytes_ou = sizeof(int)*len_wc;

      // lossless compression
      zlib_in_place_compress(qwc, &bytes_ou);

      dat_in[n_ou++] = min;
      dat_in[n_ou++] = max;
      dat_in[n_ou++] = bytes_ou;

      int_ou = (int) ceil(bytes_ou/sizeof(int))+1 ;

      // copy the data from qwc[] back into dat_in[]
      for (i=0;i<int_ou;i++)
      {
        dat_in[n_ou++]=(float) qwc[i];
      }

      wrds = wrds+3+int_ou;
   }

    par_int[5]=wrds;
    
    bytes_in = sizeof(float)*nr*nc*ns;
    cmp_factor = bytes_in/(sizeof(int)*wrds);
    printf("bytes_in= %d\t bytes_ou= %d\t cmp_factor= %f\n",bytes_in,4*wrds,cmp_factor);

    sez_free_vector(S->Win,0,(S->nc_ext)*(S->nr_ext)-1);
    sez_free_vector(S->Wou,0,(S->nc_ext)*(S->nr_ext)-1);
    sez_free_vector(S->qLA,0,S->nq-1);
    sez_free_vector(S->qLS,0,S->nq-1);
    sez_free_vector(S->qHA,0,S->nq-1);
    sez_free_vector(S->qHS,0,S->nq-1);
    sez_free_vector(qwc,0,len_wc-1); 
    free(S);
}
