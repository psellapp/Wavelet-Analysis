function [waveArray]=create_wave_array(n,dt)
% 
% Create the wavenumber (or Fourier frequencey) array for wavelet analysis. 
% The first N/2 frequencies are n=0 -> n/2; the remainder are from 
% -n/2+1 -> -1, in that order - there are a total N independent fequencies.
% 
% created: Prabu, 8/12/2015. Based on GRS' version of wv1d_e.f 
% 

waveArray = 1:fix(n/2);
waveArray = waveArray.*(2*pi()/(n*dt));
waveArray = [0, waveArray, -waveArray(fix((n-1)/2):-1:1)];
end
