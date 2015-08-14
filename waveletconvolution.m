function [transformSignal] = waveletconvolution(wv_coeff,origSignal)
% 
% Perform convolution of complex signal and wavelet. Since convolution in
% time domain is multiplication in frequency domain, the code performs
% element-wise multiplication. Make sure both inputs are in frequency space
% 
% Created: Prabu, 8/13/2015
% 

transformSignal = origSignal.*conj(wv_coeff);

end