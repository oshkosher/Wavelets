function w = ucomp(w,S)
NX=S.NX;NY=S.NY;NZ=S.NZ;
LX=S.LEV_X;LY=S.LEV_Y;LZ=S.LEV_Z;
QA=S.Q_A;QS=S.Q_S;
% Unprocess z-columns
  w = iwt_1d(w,QA,QS,LZ);
  w = reshape(w,[NZ*NX,NY]);
  w = w';
% Unprocess y-columns
  w = iwt_1d(w,QA,QS,LY);
  w = reshape(w,[NY*NZ,NX]);
  w = w';
% Unprocess x-columns
  w = iwt_1d(w,QA,QS,LX);
  w = reshape(w,[NX,NY,NZ]);'iwt: all x y z columns processed'
% Copyright (c) 2014. S.E.Zarantonello