function w=comp(in_data,S) 
NX=S.NX;NY=S.NY;NZ=S.NZ;
LX=S.LEV_X;LY=S.LEV_Y;LZ=S.LEV_Z;
QA=S.Q_A;QS=S.Q_S;
% Process x-columns
  w = in_data;
  w = reshape(w,[NX,NY*NZ]);
  w = fwt_1d(w,QA,QS,LX);
  w = w';
% Process y-columns
  w = reshape(w,[NY,NZ*NX]);
  w = fwt_1d(w,QA,QS,LY);
  w = w';
% Process z-columns
  w = reshape(w,[NZ,NX*NY]);
  w = fwt_1d(w,QA,QS,LZ);'fwt: all x y z columns processed'
% Copyright (c) 2014. S.E.Zarantonello
% Comments ?  e-mail sergio@rithmica.com