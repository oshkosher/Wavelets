function wdat = fwt_1d(udat,qa,qs,lx)

%   Description :
%    Multi 1-D forward symmetric wavelet transform using 
%    symmetric filters of odd lengths.
%   Input :
%    udat = 2-D array of unprocessed data
%    nx   = x-dimension of extended data
%    qa   = analysis low pass filter (row)
%    qs   = synthesis low pass filter (row)
%    lx   = processing level 
%   Output :
%    wdat = 2-D array of transformed data
[nx,ny] = size(udat); 

wdat = zeros(nx,ny);
na = size(qa,2); ns = size(qs,2); 
% Position qs and qa filters
if na > ns
    d  = (na-ns)/2;
    qs = [zeros(1,d),qs,zeros(1,d)];
    nf = na;
  elseif ns > na
    d  = (ns-na)/2;
    qa = [zeros(1,d),qa,zeros(1,d)]
    nf = ns;
  elseif ns == na
    nf = na;
end
rd = round((nf-1)/2); mid = rd+1;               % filter "radius"
qL = qa; tmp = qs(mid+1:nf).*(-1).^[2:rd+1];    % low-pass filter  
qH = [fliplr(tmp),-qs(mid),tmp];                % high-pass filter
% Prepare for handling borders
for r=1:2:rd
    B1(r,:) = [zeros(1,r-1),qL,zeros(1,rd-r)];
    B2(r,:) = [zeros(1,r-1),qH,zeros(1,rd-r)];
end
for r=2:2:rd
    B1(r,:) = [zeros(1,r-1),qH,zeros(1,rd-r)];
    B2(r,:) = [zeros(1,r-1),qL,zeros(1,rd-r)];
end
% Fold-over for handling borders 
V1 = fliplr(B1(:,1:rd)); V2 = fliplr(B2(:,1:rd));
F1 = B1(:,rd+1:nf+rd-1)+[zeros(rd,1),V1,zeros(rd,rd-1)];
F2 = B2(:,rd+1:nf+rd-1)+[zeros(rd,1),V2,zeros(rd,rd-1)];
F2 = flipud(fliplr(F2));
% Process x-columns
for lev = 1:lx                                      % recursive levels 
    for j=1:rd/2                                    % first & last set of coefficients
     wdat(     j,:) = F1(-1+2*j,:)*udat(1:nf-1,:);     
     wdat(nx/2+j,:) = F1(   2*j,:)*udat(1:nf-1,:); 
     k = nx/2-rd/2+j;
     wdat(     k,:) = F2(-1+2*j,:)*udat(nx-nf+2:nx,:);     
     wdat(nx/2+k,:) = F2(   2*j,:)*udat(nx-nf+2:nx,:);
    end
    k = 1;
    for j=rd/2+1:nx/2-rd/2                          % main set of coefficients
      wdat(     j,:) = qL*udat(k:k+nf-1,:);
      wdat(nx/2+j,:) = qH*udat(k+1:k+nf,:);
      k = k+2;
    end  
    nx = nx/2; 
    udat = wdat(1:nx,:);
end
% Copyright (c) 2014. S.E.Zarantonello