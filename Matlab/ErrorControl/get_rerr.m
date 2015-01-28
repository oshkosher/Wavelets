function rerr=get_rerr(wo,w,nrm,ERRTYPE)
%w=uint8(wr);wo=uint8(wo);%only when original data is uint8
if ERRTYPE=='MSE'
    err=(sum((w(:)-wo(:)).^2));
elseif ERRTYPE==' L2'
    err=sqrt(sum((w(:)-wo(:)).^2));
elseif ERRTYPE==' L1'
    err=sum(abs((w(:)-wo(:))));
end
rerr=double(err)/nrm;
% S.E.Zarantonello 2015