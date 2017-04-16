%http://www.mathworks.com/help/wavelet/ref/dwt2.html
load woman;
wname = 'haar'; %different wavelets implemented info :http://www.mathworks.com/help/wavelet/ref/wfilters.html
[CA,CH,CV,CD] = dwt2(X,wname,'mode','per'); %'mode',MODE computes the wavelet decomposition with the extension mode MODE that you specify.
%Display the vertical detail image and the lowpass approximation.
subplot(211)
imagesc(CV); title('Vertical Detail Image');
colormap gray;
subplot(212)
imagesc(CA); title('Lowpass Approximation');