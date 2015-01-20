function S = s_init(BITS,CUT,NX,NY,NZ,WAVELET)

%    Hardwired parameters       
       LEV_X =   7;%-1; 
       LEV_Y =   6;%-1;
       LEV_Z =   4;%-1;
              
%    Derived parameters       
        [Q_A,Q_S]  =  qmfilter(WAVELET);             
        L_A = size(Q_A,2); L_S = size(Q_S,2);        
        L_Q = ceil(max(L_A,L_S));
        
%    Check filter lengths
        if rem(L_A,2) == 0
          error('fwt_3d error: analysis filter has even length')
          return
        elseif rem(L_S,2) == 0
          error('fwt_3d error: synthesis filter has even length')
          return
        end 
        
%      Default values of LEV_X, LEV_Y, LEV_Z
       if LEV_X == -1                                   
        LEV_X = floor(log2(NX/L_Q))+1;
       end       
       if LEV_Y == -1                                  
        LEV_Y = floor(log2(NY/L_Q))+1; 
       end  
       if LEV_Z == -1                                  
        LEV_Z = floor(log2(NZ/L_Q))+1; 
       end   
                   
%   Build structure
       S = struct('BITS',BITS,'CUT',CUT,'LEV_X',LEV_X,'LEV_Y',LEV_Y,...
                  'LEV_Z',LEV_Z,'NX',NX,'NY',NY,'NZ',NZ,'Q_A',Q_A,'Q_S',Q_S);  
% Copyright (c) 2014. S.E.Zarantonello