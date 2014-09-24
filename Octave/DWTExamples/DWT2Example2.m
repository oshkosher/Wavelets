%http://stackoverflow.com/questions/6766210/matlab-wavelet-transform-for-n-dimensionmodified

X = imread('peppers.png');  %# Load sample image
nColors = 256;              %# Number of values per color component
nColors = double(intmax(class(X)))+1;
nLevel = 3;             %# Number of decompositions
cA = cell(1,nLevel);    %# Approximation coefficient storage
cH = cell(1,nLevel);    %# Horizontal detail coefficient storage
cV = cell(1,nLevel);    %# Vertical detail coefficient storage
cD = cell(1,nLevel);    %# Diagonal detail coefficient storage
startImage = X;
for iLevel = 1:nLevel,  %# Apply nLevel decompositions
  [cA{iLevel},cH{iLevel},cV{iLevel},cD{iLevel}] = dwt2(startImage,'db1');
  startImage = cA{iLevel};
end

%This is just to display the image
tiledImage = wcodemat(cA{nLevel},nColors);
for iLevel = nLevel:-1:1
  tiledImage = cat(1,cat(2,tiledImage,...
                           wcodemat(cH{iLevel},nColors)),...
                     cat(2,wcodemat(cV{iLevel},nColors),...
                           wcodemat(cD{iLevel},nColors)));
end
figure;
imshow(uint8(tiledImage-1));  %# Convert to unsigned 8-bit integer to display