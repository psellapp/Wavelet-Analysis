% script to perform wavelet analysis on signal in variable y. Load/create
% variable y prior to running this script
% 
% Created: Prabu, 8/13/2015
% 

linlog_flag = 0; % 0 - linear scale; 1 - log scale
origSignal = fft(y);
n = length(origSignal);

w0 = 7.5; %Morlet wavelet parameter. 

[waveArray]=create_wave_array(n,[]);
[scale] = waveletscale(n,linlog_flag);
transformSignal_array = zeros(length(scale),n);

for i = 1:length(scale)
    [Morlet_hat_coeff]=Morlet_hat(waveArray,scale(i),w0);
    [transformSignal] = waveletconvolution(Morlet_hat_coeff,origSignal);
    transformSignal_array(i,:) = ifft(transformSignal);
end

contourf(t,scale,abs(transformSignal_array))