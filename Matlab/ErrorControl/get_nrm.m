function nrm=get_nrm(w,ERRTYPE) 
if ERRTYPE=='MSE'
    nrm=(sum(w(:).^2));
elseif ERRTYPE==' L2'
    nrm=sqrt(sum(w(:).^2));
elseif ERRTYPE==' L1'
    nrm=sum(abs(w(:)));
end
% S.E.Zarantonello 2015