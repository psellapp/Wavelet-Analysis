function [transformSignal_array] = wavelet_run_script(test_signal,x_var,min_scale,max_scale,ds)
% Function to perform wavelet analysis
%
% Created: Prabu, 8/13/2015
% modified: Prabu, 8/24/2015. -generate test signal
%           Prabu, 9/3/2015. -calculate energy
%           Prabu, 4/20/2016. -Major changes, code actually works now! Call as a function, outputs in terms of scales themselves  

%================= Generate test signal====================================
%============= Comment out appropriately to generate req signal============

% Fs = 8000; %Sampling rate, Hz
% dt = 1/Fs;
% t = 0:dt:1-dt;
% f1 = 100;
% f2 = 34; %test frequencies, in Hz
% 
% test_signal1 = sin(2*pi()*f1*t);
% test_signal2 = 0.*t;
% 
% test_signal1 = sin(2*pi()*f1*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3)+4*rand(1,length(t)).*(t>0.7);
% test_signal2 = 0.*t;

% test_signal1 = (sin(2*pi()*f1*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3))+rand(1,length(t));
% test_signal2 = (sin(2*pi()*(f1+10)*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3));
% 
% test_signal1 = sin(2*pi()*f1*t).*(t<=0.3)+sin(2*pi()*f2*t).*(t>0.3);
% test_signal2 = (sin((2*pi()*f1*t)+(pi/2)).*(t<=0.3))+(sin(2*pi()*f2*t).*(t>0.3));

% test_signal1 = sin(2*pi()*f1*t);
% test_signal2 = sin((2*pi()*f2*t));
% 

% test_signal = sin(2*pi()*4*t).*(t<0)+sin(2*pi()*1*t).*(t>=0);
%  test_signal = test_signal1 + test_signal2;
% % test_signal = test_signal .*1;
% test_signal = rho1_Nsquared_filterFc300(1:end)';
% clear t;
%  t = linspace(0,1,length(test_signal));
%  dt = t(2)-t(1);
%==========================================================================
% ====== calculate params for real signal =================================
% test_signal = (rho5_fit_filt300 - rho5_sorted)';
% R = 4; %radius of grid in cm
% t = z_norm.*R; % height in cm

% min_scale = 0.05; % specify in cm
% max_scale = 1.5; %specify in cm

t = x_var;
test_signal = test_signal';

dt = abs(t(2) - t(1)); % height increment 
% ds = .005; %scale increment - make sure it is greater than dt
linlog_flag = 0; % 0 - linear scale; 1 - log scale
% =========================================================================

test_signal_fft = fft(test_signal);
n = length(test_signal_fft);
signal_energy = sum(test_signal.^2)/n;%calculate energy in test signal
e_fft = sum(abs(test_signal_fft).^2)/n^2;
% =========================================================================



[waveArray]=create_wave_array(n,dt);



a0 = min_scale;
a1 = max_scale;

[scale,delta,da] = waveletscale(n,dt,ds,a0,a1,linlog_flag);
%  scale = 2*pi/n:.001:pi;
transformSignal_array = zeros(length(scale),n);
cg = 0;
% enorm = zeros(1,n);
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
%         morlet_fourier_factor = 4*pi()/(w0+sqrt(2+w0^2));
%         freqScale = 1./(morlet_fourier_factor.*scale);
%         pseudo_freq = freqScale;
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
figure(5)
plot(scale,e_wv,'+r')
xlabel('Scale')
ylabel('Energy')
%==========================================================================
% figure(15)
% xScale = t;
% subplot(3,1,3)
% contourf((xScale),scale,abs(transformSignal_array))
% ylabel('Scale, cm');


figure(9)
xScale = t;

subplot(3,1,1)
plot((xScale),test_signal)
ylabel('\rho '', gm/cm^3')
axis tight

subplot(3,1,2)
% contourf(xScale,pseudo_freq,angle(transformSignal_array).*(abs(transformSignal_array)>0.005))
% contourf(xScale,pseudo_freq,abs(transformSignal_array))
% contourf(fliplr(xScale),scale./w0,abs(transformSignal_array))
contourf((xScale),scale,abs(transformSignal_array))
ylabel('Scale, cm');

subplot(3,1,3)
% contourf(fliplr(xScale),scale./w0,angle(transformSignal_array).*(abs(transformSignal_array) > 1.5e-6))
contourf((xScale),scale,angle(transformSignal_array).*(abs(transformSignal_array) > 9e-6))

ylabel('Scale, cm');
xlabel('Z, cm')
% contourf(xScale,pseudo_freq,abs(transformSignal_array))
% contourf(xScale,pseudo_freq,abs(transformSignal_array).*(abs(transformSignal_array)>0.3))
% 
% figure(15)
% [reconstructedSignal] = reconstructSignal(transformSignal_array,scale,waveArray,w0);
% plot(reconstructedSignal,'b-')
% hold on
% plot(test_signal,'r--')
% hold off