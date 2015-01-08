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
    int     nr, nc;  // number of rows, number of columns
    int     nq;      // filter length (currently fixed at 9)

    int     nr_ext, nc_ext;  // number of row,columns after expansion

    int     ir_0, ic_0;  // index of the first element in each row, col
    int     ir_1, ic_1;  // 0=expanded data, 1=initial data
    int     jr_0, jc_0;  // index of the last element in each row, col
    int     jr_1, jc_1;

    int     lev_r, lev_c;  // # of transform levels for rows, columns

    float*  qLA;   // arrays of filter constants for wavelet convolution
    float*  qHA;   // See "Hard-wired analysis and synthesis filters"
    float*  qLS;   // in sez_wcomp.c
    float*  qHS;

    float*  Win;
    float*  Wou;
} sez_wcomp_par;

sez_wcomp_par* sez_wcomp_init();
void sym_fwt_2d();
void sym_iwt_2d();
#endif /*__SEZ_WCOMP_H__*/
