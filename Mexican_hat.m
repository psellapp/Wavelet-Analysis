function [Mexican_hat_coeff]=Mexican_hat(waveArray,scale)
%
% Calculate coefficients for the transform of Mexican Hat wavelet.
% The results are returned as complex numbers even for real-valued
% functions.
%
% created: Prabu, 8/16/2015. Based on GRS' version of Mexican_hat
%
% waveArray - 1-D vector of wavenumbers. Generate using function create_wave_array

rp = zeros(1,length(waveArray));
rp(waveArray>=0) = exp(-((scale*waveArray(waveArray>=0)).^2)/2).*((scale*waveArray(waveArray>=0)).^2)*(-2);

Mexican_hat_coeff = complex(rp);

% plot(abs(Mexican_hat_coeff));

end