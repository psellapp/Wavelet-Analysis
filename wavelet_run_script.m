% script to perform wavelet analysis on test signal 
%
% Created: Prabu, 8/13/2015
% modified: Prabu, 8/24/2015. -generate test signal

clear
%================= Generate test signal====================================
Fs = 1000; %Sampling rate, Hz
dt = 1/Fs;
t = 0:dt:1-dt;
f1 = 20; 
f2 = 30; %test frequencies, in Hz
test_signal = sin(2*pi()*f1*t);%.*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3);
test_signal_fft = fft(test_signal);
% =========================================================================

linlog_flag = 0; % 0 - linear scale; 1 - log scale
n = length(test_signal_fft);

[waveArray]=create_wave_array(n,[]);

a0 = [];
a1 = [];
ds = .01; 
[scale] = waveletscale(n,dt,ds,a0,a1,linlog_flag);
%  scale = 2*pi/n:.001:pi;
transformSignal_array = zeros(length(scale),n);


% =========================================================================
% Select wavelet - comment/uncomment as required

% =============== Morlet wavelet ==========================================
wavelet_name = 'morlet';
w0 = 5.5; %Morlet wavelet parameter.

% =============== Mexican Hat wavelet =====================================
% wavelet_name = 'mexican_hat';
% =========================================================================

if strcmp(wavelet_name,'mexican_hat')
    disp('Wavelet selected - Mexican Hat')
    for i = 1:length(scale)
        [Mexican_hat_coeff]=Mexican_hat(waveArray,scale(i));
        [transformSignal] = waveletconvolution(Mexican_hat_coeff,test_signal_fft);
        transformSignal_array(i,:) = ifft(transformSignal);
    end
    pseudo_freq = scale;
elseif strcmp(wavelet_name,'morlet')
    disp('Wavelet selected - Morlet')
    for i = 1:length(scale)
        [Morlet_hat_coeff]=Morlet_hat(waveArray,scale(i),w0);
        [transformSignal] = waveletconvolution(Morlet_hat_coeff,test_signal_fft);
        transformSignal_array(i,:) = ifft(transformSignal);
    end
%     morlet_fourier_factor = 4*pi()/(w0+sqrt(2+w0^2));
%     freqScale = 1./(morlet_fourier_factor.*scale); 
    pseudo_freq = w0./scale;%scale to pseudo-freq. w0 is center freq.
else
    disp('Exciting new wavelets coming Fall 2015!')
end

figure(1)
xScale = t;
subplot(2,1,1)
plot(xScale,test_signal)
ylabel('Original signal, y, [amplitude]')
subplot(2,1,2)
contourf(xScale,pseudo_freq,abs(transformSignal_array))
ylabel('Pseudo-frequency, [Hz]');
xlabel('Time [s]')