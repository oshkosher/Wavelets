        #include <stdio.h>
        #include <math.h>
	#include "sez_alloc.h" 

        //../TEST12/SourceWave.wi.le.H@
        //out_n.dat
	#define nz 16 // depth
	#define nx 16 // xline 
	#define ny  1 // inline 
	#define nt    1 // time steps

        /*
        //../TEST12/surfzx.H@
	#define nz  100 // depth
	#define nx  100 // inline 
	#define ny    1 // xline 
	#define nt    1 // time steps

        //../TEST12/TT.H@
	#define nz  201 // depth
	#define nx  101 // inline 
	#define ny    1 // xline 
	#define nt    1 // time steps
        */
       
void comp(float*,int*);

int	main()
	{
	int  i,ix,iz,it,n,qlen,wrds,*par_int,*wrd_int;
        float cmp_ratio,*dat_in,*par_float;
        FILE *fptr,*fcmp;
        
/******* CREATE INPUT DATA BUFFER ******/
        dat_in    = sez_float_vector(0,nz*nx*ny-1); //allocate 1-D array

        /* 
         for (i=0;i<nx*ny*nz;i++)
         {
          dat_in[i] = 0.;
         }  
        */ 

/******* OPEN OUTPUT COMPRESSED DATA FILE ******/
        if (( fcmp = fopen("cmp.dat","w")) == NULL){
          printf(" ERROR IN COMP FILE OPEN !!! ");
          return 0;
        } 

/******* SAVE SPACE FOR PREAMBLE IN OUTPUT FILE ******/
        par_int   = sez_int_vector(0,5); //allocate parameter array
        wrd_int  = sez_int_vector(0,nt-1); //allocate "size" array

        fwrite(par_int,4,5,fcmp);  //will overwrite later
        fwrite(wrd_int,4,nt,fcmp); //will overwrite later

/******* LOAD PARAMETERS ******/
        par_int[0] = nz; 
        par_int[1] = nx; 
        par_int[2] = ny; 
        par_int[3] = 12000; //comp. par.(bigger=less err & less compress)
        par_int[4] = nt;   //number of snapshots 
        par_int[5] = 0;   //returns the number of words in comp. buffer 

/******* TIME-SNAPSHOT LOOP ******/

        for (it=0;it<nt;it++)
        {

        //wavefield buffer updated here(not done in this example for simplicity) 
        //if (( fptr = fopen("out_t.dat","r")) == NULL){
        //if (( fptr = fopen("../TEST12/TT.H@","rl")) == NULL){
        //if (( fptr = fopen("../TEST12/surfzx.H@","rl")) == NULL){
        // if (( fptr = fopen("../TEST12/SourceWave.wi.le.H@","rl")) == NULL){
        if (( fptr = fopen("sez.in","rb")) == NULL){
          printf(" ERROR IN FILE OPEN !!! ");
          return 0;
        } else
	{
          fread(dat_in,4,nz*nx*ny,fptr);
        } 
        fclose(fptr);      
        
        /** COMPRESS WAVEFIELD SNAPSHOT **/
           comp(dat_in,par_int);

           wrds = par_int[5]; //number of words in "snapshot" compressed buffer

           wrd_int[it] = wrds; 

        /** WRITE "SNAPSHOT" COMPRESSED BUFFER TO DISK **/
           fwrite(dat_in,4,wrds,fcmp);
        }

/******* REPOSITION & WRITE PREAMBLE TO DISK ******/
        fseek(fcmp,0,0);      
        fwrite(par_int,4,5,fcmp); 
        fwrite(wrd_int,4,nt,fcmp); 

/******* FREE MEMORY & CLOSE COMPRESSED DATA FILE ******/
        fclose(fcmp);      
        sez_free_vector(dat_in,0, nz*nx*ny-1);
        sez_free_int(par_int,0,5);
        sez_free_int(wrd_int,0,nt-1);

        }
