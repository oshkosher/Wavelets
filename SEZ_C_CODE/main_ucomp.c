        #include <stdio.h>
        #include <math.h>
	#include "sez_alloc.h" 

void ucomp(float*,float*,int*);

int	main()
	{
	int i,it,nx,ny,nz,nt,wrds,start,*par_int,*wrd_int;
        float *dat_re,*dat_ou;
        FILE *fcmp,*fout;
        
/******* ALLOCATE MEMORY ******/
        par_int   = sez_int_vector(0,5);

/******* OPEN COMPRESSED DATA FILE ******/
        if (( fcmp = fopen("cmp.dat","r")) == NULL){
        //if (( fcmp = fopen("SourceWave","r")) == NULL){
          printf(" ERROR IN COMP FILE OPEN !!! ");
          return 0;
        } 

/******* SPECIFY TIME SLICE ******/
        it = 0;

/******* READ PARAMETERS ******/
        fread(par_int,4,5,fcmp);
        nz   = par_int[0]; 
        nx   = par_int[1];
        ny   = par_int[2];
        nt   = par_int[4];

        wrd_int  = sez_int_vector(0,nt-1);

        fread(wrd_int,4,nt,fcmp);

        wrds = wrd_int[it];
        par_int[5] = wrds;
        start= sizeof(int)*(5+nt); // bytes in prologue        
        for (i=0;i<it;i++){start+=sizeof(int)*wrd_int[i];} 
 
/******* ALLOCATE ARRAYS, READ COMPRESSED FILE & UNCOMPRESS ******/
        dat_re  = sez_float_vector(0,nz*nx*ny-1);
        dat_ou  = sez_float_vector(0,wrds-1);

        //fseek(fcmp,4*(wrd_int[0]+nt+5),0);
        fseek(fcmp,start,0);
        fread(dat_ou,4,wrds,fcmp);
        ucomp(dat_ou,dat_re,par_int);
       
//******* WRITE RECREATED DATA ******/
        if (( fout = fopen("out.dat","w")) == NULL){
          printf(" ERROR IN OUTPUT FILE OPEN !!! ");
          return;
        }
        fwrite(dat_re,4,nx*ny*nz,fout);

/******* FREE MEMORY & CLOSE FILES ******/
        sez_free_vector(dat_ou,0, wrds-1);
        sez_free_vector(dat_re,0, nz*nx*ny-1);

        fclose(fcmp);
        fclose(fout);
        }
