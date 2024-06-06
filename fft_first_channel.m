% Compute FFT of the first channel
% data =  load('FN_032_V1_Postdose1_Nback.mat');
% intensityData = data.d;
% 
% firstChannel = intensityData(:, 1);
n = length(firstChannel);
% Sampling frequency
fs = 1 / (timeVector(2) - timeVector(1)); 
% Normalize FFT by the sample size
y = fft(firstChannel) / n; 
% Frequency vector
f = (0:n-1) * (fs / n); 

% Plot the Fourier transform
figure;
plot(f, abs(y)*2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0.5 4])
title('Fourier Transform of the First Channel');

% Approximate pulse frequency range in Hz (1-2 Hz)
pulseFreqRange = [1, 2]; 
% Frequencies from 2.5 Hz to the max for the noise calculation as asked in
% the hw problem
noiseFreqRange = [2.5, max(f)]; 

% Find the index of the pulse frequency range
pulseFreqIdx = (f >= pulseFreqRange(1)) & (f <= pulseFreqRange(2));
% Find the peak within the pulse frequency range
[~, signalIdx] = max(abs(y(pulseFreqIdx))); 
% Find that peak as is in the ENTIRE signal
signalIdx = find(pulseFreqIdx, 1) + signalIdx - 1;

% Calculate the noise level
noiseFreqIdx = (f >= noiseFreqRange(1)) & (f <= noiseFreqRange(2));
% Average noise level in the specified range of noise
noise = mean(abs(y(noiseFreqIdx)));
% Signal strength at the pulse frequency - which is the peak within the
% range of the heartbeat
signal = abs(y(signalIdx)); 

% Compute SNR - signal to noise ratio
SNR = signal / noise;
disp(SNR);
