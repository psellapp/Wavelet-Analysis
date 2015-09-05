% script to perform wavelet analysis on test signal
%
% Created: Prabu, 8/13/2015
% modified: Prabu, 8/24/2015. -generate test signal
%           Prabu, 9/3/2015. -calculate energy

clear
%================= Generate test signal====================================
%============= Comment out appropriately to generate req signal============

Fs = 1024; %Sampling rate, Hz
dt = 1/Fs;
t = 0:dt:1-dt;
f1 = 20;
f2 = 45; %test frequencies, in Hz

test_signal1 = sin(2*pi()*f1*t);
test_signal2 = 0.*t;

% test_signal1 = (sin(2*pi()*f1*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3))+rand(1,length(t));
% test_signal2 = (sin(2*pi()*(f1+10)*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3));

% test_signal1 = sin(2*pi()*f1*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3);
% test_signal2 = (sin((2*pi()*f1*t)+(pi/2)).*(t<=0.3))+(sin(2*pi()*f2*t).*(t>0.3));

% test_signal1 = sin(2*pi()*f1*t);
% test_signal2 = sin((2*pi()*f2*t)+(pi/2));

test_signal = test_signal1 + test_signal2;
%==========================================================================


test_signal_fft = fft(test_signal);
n = length(test_signal_fft);
signal_energy = sum(test_signal.^2)/n;%calculate energy in test signal
e_fft = sum(abs(test_signal_fft).^2)/n^2;
% =========================================================================

linlog_flag = 1; % 0 - linear scale; 1 - log scale

[waveArray]=create_wave_array(n,[]);

a0 = [10];
a1 = [];
ds = .01;
[scale,delta,da] = waveletscale(n,dt,ds,a0,a1,linlog_flag);
%  scale = 2*pi/n:.001:pi;
transformSignal_array = zeros(length(scale),n);
cg = 0;
enorm = zeros(1,n);
emod = zeros(1,n);
e_wv = zeros(1,length(scale));

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
        
        cg = cg + sum(delta.*abs(Mexican_hat_coeff).^2);%calculate integration constant for given wavelet function. Admissibility condition
        enorm = delta.*(abs(transformSignal_array(i,:)).^2);
        emod = emod + enorm;
        e_wv(i) = sum(enorm);
        
    end
    pseudo_freq = scale;
elseif strcmp(wavelet_name,'morlet')
    disp('Wavelet selected - Morlet')
    for i = 1:length(scale)
        [Morlet_hat_coeff]=Morlet_hat(waveArray,scale(i),w0);
        [transformSignal] = waveletconvolution(Morlet_hat_coeff,test_signal_fft);
        transformSignal_array(i,:) = ifft(transformSignal);
        
        cg = cg + sum(delta.*abs(Morlet_hat_coeff).^2);%calculate integration constant for given wavelet function. Admissibility condition
        enorm = delta.*(abs(transformSignal_array(i,:)).^2);
        emod = emod + enorm;
        e_wv(i) = sum(enorm);
    end
    %     morlet_fourier_factor = 4*pi()/(w0+sqrt(2+w0^2));
    %     freqScale = 1./(morlet_fourier_factor.*scale);
    pseudo_freq = w0./scale;%scale to pseudo-freq. w0 is center freq.
else
    disp('Exciting new wavelets coming Fall 2015!')
end

%=============Normalize energy quantities by cg and print to screen========
e_wv = e_wv./cg;%Energy at each scale
etot = sum(emod)/cg;
disp('Energy in wavelet. Cg = ')
cg
disp('Total energy after wavelet transform:')
etot
disp('Energy in original signal')
signal_energy
disp('Energy from PSD of original signal')
e_fft
figure(4)
plot(scale,e_wv,'+r')
xlabel('Scale')
ylabel('Energy')
%==========================================================================

figure(1)
xScale = t;
subplot(2,1,1)
plot(xScale,test_signal)
ylabel('Original signal, y, [amplitude]')
subplot(2,1,2)
contourf(xScale,pseudo_freq,angle(transformSignal_array).*(abs(transformSignal_array)>0.05))
ylabel('Pseudo-frequency, [Hz]');
xlabel('Time [s]')