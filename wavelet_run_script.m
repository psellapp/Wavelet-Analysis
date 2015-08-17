% script to perform wavelet analysis on signal in variable y. Load/create
% variable y prior to running this script. Variable should contain 1-D data
%
% Created: Prabu, 8/13/2015
%

linlog_flag = 0; % 0 - linear scale; 1 - log scale
origSignal = fft(y);
n = length(origSignal);
[waveArray]=create_wave_array(n,[]);
[scale] = waveletscale(n,linlog_flag);
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
        [transformSignal] = waveletconvolution(Mexican_hat_coeff,origSignal);
        transformSignal_array(i,:) = ifft(transformSignal);
    end
elseif strcmp(wavelet_name,'morlet')
    disp('Wavelet selected - Morlet')
    for i = 1:length(scale)
        [Morlet_hat_coeff]=Morlet_hat(waveArray,scale(i),w0);
        [transformSignal] = waveletconvolution(Morlet_hat_coeff,origSignal);
        transformSignal_array(i,:) = ifft(transformSignal);
    end
else
    disp('Exciting new wavelets coming Fall 2015!')
end

xScale = linspace(0,1,length(y));
subplot(2,1,1)
plot(xScale,y)
ylabel('Original signal, y')
subplot(2,1,2)
contourf(xScale,scale,abs(transformSignal_array))
ylabel('Scale');