/**
 * \file sez_wcomp.h
 * 
 * TT lossy compression
 *
 *
 * $Date: 2006/06/13 17:00:00 $
 *
 * \author Sergio E. Zarantonello
 */
#ifndef __SEZ_WCOMP_H__     
#define __SEZ_WCOMP_H__     

typedef struct sez_wcomp_par_s {
    int     nq, nr, nc, nr_ext, nc_ext;  
    int     ir_0, ic_0;
    int     ir_1, ic_1;
    int     jr_0, jc_0;
    int     jr_1, jc_1;
    int     lev_r, lev_c;
    float*  qLA;
    float*  qHA;
    float*  qLS;
    float*  qHS;
    float*  Win;
    float*  Wou;
} sez_wcomp_par;

sez_wcomp_par* sez_wcomp_init();
void sym_fwt_2d();
void sym_iwt_2d();
#endif /*__SEZ_WCOMP_H__*/
