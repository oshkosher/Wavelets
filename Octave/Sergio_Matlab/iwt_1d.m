function udat = iwt_1d(wdat,qa,qs,lx)

%   Description :
%    multiple 1-D inverse symmetric wavelet transform using 
%    symmetric filters of odd length.
%   Input :
%    wdat = 2-D array of processed data
%    qa   = analysis low pass filter (row)
%    qs   = synthesis low pass filter (row)
%    lx   = processing level in x-direction
%   Output :
%    udat = 2-D inversed data

[nx,ny] = size(wdat); udat = wdat;
mx=nx;

na = length(qa); ns = length(qs);

% Create inverse filters
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

rd  = round((nf-1)/2); mid = rd+1;                               % filter "radius"
qL  = qs; tmp = qa(mid+1:nf).*(-1).^[2:rd+1];                    % low-pass filter  
qH  = [fliplr(tmp),-qa(mid),tmp];                                % high-pass filter
tLH = [qL(1:2:nf);[qH(2:2:nf),0]]; tLH=tLH(:); qLH=tLH(1:nf);    % interleaved LH filter
tHL = [qH(1:2:nf);[qL(2:2:nf),0]]; tHL=tHL(:); qHL=tHL(1:nf);    % interleaved HL filter

for r=1:2:rd
    B1(:,r) = [zeros(r-1,1);qLH;zeros(rd-r,1)];
    B2(:,r) = [zeros(r-1,1);qHL;zeros(rd-r,1)];
end
for r=2:2:rd
    B1(:,r) = [zeros(r-1,1);qHL;zeros(rd-r,1)];
    B2(:,r) = [zeros(r-1,1);qLH;zeros(rd-r,1)];
end

% Fold-over

V1 = flipud(B1(1:rd,:)); V2 = flipud(B2(1:rd,:));
F1 = B1(rd+1:nf+rd-1,:)+[zeros(1,rd);V1;zeros(rd-1,rd)];
F2 = B2(rd+1:nf+rd-1,:)+[zeros(1,rd);V2;zeros(rd-1,rd)];
F2 = fliplr(flipud(F2));
F1 = F1'; F2 = F2';

% Process x-columns

for lev = 1:lx                                            % recursive levels 
    lnx = 2^(lx-lev); nx = mx/lnx;
    wdat(1:nx,:) = [udat(1:nx/2,:); wdat(nx/2+1:nx,:)];
    reorder = [1:nx/2;nx/2+1:nx]; reorder = reorder(:);   % reorder data
    tmp = wdat(reorder,:);
    for j=1:rd                                            % first & last batch of coefficients
     udat(j,:) = F1(j,:)*tmp(1:nf-1,:);     
     k = nx-rd+j;
     udat(k,:) = F2(j,:)*tmp(nx-nf+2:nx,:);     
    end
    k = 1;
    for j=rd+1:2:nx-rd                                    % main body of coefficients
      udat(  j,:) = qLH'*tmp(k:k+nf-1,:);
      udat(1+j,:) = qHL'*tmp(k+1:k+nf,:);
      k = k+2;
    end   
end
wdat = udat;
% Copyright (c) 2013. S.E.Zarantonello