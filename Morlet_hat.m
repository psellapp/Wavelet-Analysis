function [Morlet_hat_coeff]=Morlet_hat(waveArray,scale,w0)
%
% Calculate coefficients for the non-analytic transform of Morlet wavelet.
% The results are returned as complex numbers even for real-valued
% functions.
%
% created: Prabu, 8/12/2015. 
% Based on GRS' version of Morlet_hat and wave_bases code from Christopher
% Torrence and Gilbert P. Compo, University of Colorado
%
% waveArray - vector of wavenumbers. Generate using function create_wave_array
% w0 - morlet hat parameter. Defaults to 5.5 if empty
% 

if isempty(w0)
    w0 = 5.5;
end

norm = sqrt(scale*waveArray(2)*length(waveArray))/pi^(0.25); %normalize wavelet so that sum(Morlet_hat_coeff.^2)/length(waveArray) = 1
arg = norm*exp(-(((scale.*waveArray)-w0).^2/2).*(waveArray > 0)); 
Morlet_hat_coeff = complex(arg.*(waveArray > 0)); 

%   figure(3);
%   plot(real(ifft(Morlet_hat_coeff)).*(waveArray>0)+real(ifft(Morlet_hat_coeff)).*(waveArray<=0));
% 
end