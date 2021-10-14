%carrier frequency:
carrierFreq = 10000; %10kHz for carrier frequency

%carrier signal 16 times oversampled:
samplingFreq = 16 * carrierFreq; %sampling frequency is 16 times the carrier frequency

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024; %defined from Phase 1

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%define amplitude of 2 for what??
amplitude = 8; 

%time scale in seconds for ....
time = bits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%Assume a 6th order BANDPASS filter with cut-off frequency 0.2:
[b, a] = butter(6, 0.2); %low-pass filter

%carrier signal
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);

%signal length is for everywhere
signalLength = samplingFreq * bits/dataRate + 1;

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of run of times
runTime = 10;

ookErrorRate = zeros(length(SNR)); %OOK error rate








