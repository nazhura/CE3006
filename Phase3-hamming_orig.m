%
% Data Assumptions
% 
% Amplitude of OOK is set to 8
% Hamming Code: Codeword 7 bits, Data 4 bits
% Sampling Frequency = 16 x Carrier Frequency
% 
% Generate data -> Hamming Code Encoding -> OOK/BPSK Modulation -> Noise -> Demodulation -> Comparison
% 

%carrier frequency:
carrierFreq = 10000; %10kHz
carrierFreqFSK1 = 10000;
carrierFreqFSK2 = 40000;

%Self-defined: Codeword length (n) & Message length (k)
codeword_length = 7;
message_length = 4;

%carrier signal 16 times oversampled:
samplingFreq = 16* carrierFreq; %sampling frequency is 16 times the carrier frequency

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;
encoded_bits = bits*codeword_length/message_length;

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%amplitude for
amplitude = 8;

%timescale in seconds for ....
time = encoded_bits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%assume a 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%Carrier Signal Generation
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);
carrierSignal_FSK1 = amplitude .* cos(2*pi*carrierFreqFSK1*timeScale);
carrierSignal_FSK2 = amplitude .* cos(2*pi*carrierFreqFSK2*timeScale);

%signal length is for everywhere
signalLength = samplingFreq*encoded_bits/dataRate + 1;
orig_SignalLength = samplingFreq* bits/dataRate + 1;

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 100;

ookErrorRate = zeros(length(SNR));
ookErrorRate_org = zeros(length(SNR));
bpskErrorRate = zeros(length(SNR));
bfskErrorRate = zeros(length(SNR));
        
%generate data - Hamming
dataHamming = round(rand(1, bits));
dataHamming = transpose(dataHamming);

%encoding - hamming
encodedHamming = encode(dataHamming, codeword_length, message_length ,'hamming/binary');
encodedHamming = transpose(encodedHamming);

%sampling
dataSignal = zeros(1, signalLength);
actualDataSignal = zeros(1, orig_SignalLength);
for k = 1: signalLength - 1
    signal(k) = encodedHamming(ceil(k*dataRate/samplingFreq));
end
for k = 1: orig_SignalLength - 1
    actualDataSignal(k) = dataSignal(ceil(k*dataRate/samplingFreq));
end
dataSignal(signalLength) = dataSignal(signalLength - 1);
actualDataSignal(orig_SignalLength) = actualDataSignal(orig_SignalLength - 1);

%==== OOK ====%
%encoded signal
ookSignal = carrierSignal .* dataSignal;
ookSignalPower = (norm(ookSignal)^2)/signalLength;
ookNoisePower = ookSignalPower ./ SNR;

%unencoded signal
orig_ookSignal = carrierSignal(1:orig_SignalLength) .* actualDataSignal;
orig_ookSignalPower = (norm(orig_ookSignal)^2)/orig_SignalLength;
orig_avgOOKNoisePower = orig_ookSignalPower ./ SNR;

%==== BPSK ====%
bpskSourceSignal = dataSignal .* 2 - 1;
bpskSignal = carrierSignal .* bpskSourceSignal;
bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
bpskNoisePower = bpskSignalPower ./ SNR;

%==== BFSK ====%
bfskSourceSignal_high = carrierSignal_FSK1 .* (dataSignal == 1);
bfskSourceSignal_low = carrierSignal_FSK2 .* (dataSignal == 0);
bfskSignal = bfskSourceSignal_low + bfskSourceSignal_high;
bfskSignalPower = (norm(bfskSignal)^2)/signalLength;
avgBFSKNoisePower = bfskSignalPower ./ SNR;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;
    
	for j = 1 : testSamples
      
        % encoded OOK
        avgOOKNoisePower = ookSignalPower ./SNR(i);
        ookNoise = sqrt(avgOOKNoisePower) .*transpose(randi([0 1], signalLength, 1));
        receivedOOKSignal = ookSignal + ookNoise;
        
        % unencoded OOK
        orig_avgOOKNoisePower = orig_ookSignalPower ./ SNR(i);
        orig_ookNoise =  sqrt(orig_avgOOKNoisePower) .*transpose(randi([0 1],orig_SignalLength, 1));
        orig_receivedOOKSignal = orig_ookSignal + orig_ookNoise;
        
        % BPSK
        avgBPSKNoisePower = bpskSignalPower ./SNR(i);
        bpskNoise = sqrt(avgBPSKNoisePower) .*transpose(randi([0 1],signalLength, 1));
        receivedBPSKSignal = bpskSignal + bpskNoise;
        
        % BFSK
        avgBFSKNoisePower = bfskSignalPower ./SNR(i);
        bfskNoise = sqrt(avgBFSKNoisePower) .*transpose(randi([0 1],signalLength, 1));
        receivedBFSKSignal = bfskSignal + bfskNoise;
        
        % Non-coherent Detection method of demodulation
        ookDemodulated = receivedOOKSignal .* 2.* carrierSignal;
        ookFiltered = filtfilt(b, a, ookDemodulated);
        
        orig_ookDemodulated = orig_receivedOOKSignal .* (2 .* carrierSignal(1:orig_SignalLength));
        orig_ookFiltered = filtfilt(b, a, orig_ookDemodulated);        
        
        bpskDemodulated = receivedBPSKSignal .* (2.* carrierSignal);
        bpskFiltered = filtfilt(b, a, bpskDemodulated);
        
        %BFSK detection
        bfskDemodulated_high = receivedBFSKSignal .* (2 .* carrierSignal_FSK1);
        bfskFiltered_high = filtfilt(b, a, bfskDemodulated_high);
        bfskDemodulated_low = receivedBFSKSignal .* (2 .* carrierSignal_FSK2);
        bfskFiltered_low = filtfilt(b, a, bfskDemodulated_low);
        bfskFiltered = bfskFiltered_high - bfskFiltered_low;
        
        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;
        
        [ookInput, ookOutput] = samplingforthreshold(ookFiltered, samplingPeriod, avgPower, encoded_bits);
        [orig_ookInput, orig_ookOutput] = samplingforthreshold(orig_ookFiltered, samplingPeriod, avgPower, bits);
        [bpskInput, bpskOutput] = samplingforthreshold(bpskFiltered, samplingPeriod, 0, encoded_bits);
        [bfskInput, bfskOutput] = samplingforthreshold(bfskFiltered, samplingPeriod, 0, encoded_bits);

        %Cyclic Decoding
        decoded_OOK = decode(ookOutput,codeword_length, message_length,'hamming/binary');
        decoded_BPSK = decode(bpskOutput,codeword_length, message_length,'hamming/binary');
        decoded_BFSK = decode(bfskOutput,codeword_length, message_length,'hamming/binary');
        
        ookError =  biterr(decoded_OOK, generatedData) ./bits;
        orig_ookError = biterr(orig_ookOutput, generatedData) ./bits;
        bpskError = biterr(decoded_BPSK, generatedData) ./bits;
        bfskError = biterr(decoded_BFSK, generatedData) ./bits;
        
        ookAvgError = ookError + ookAvgError;
        orig_ookAvgError = orig_ookError + orig_ookAvgError;
        bpskAvgError = bpskError + bpskAvgError;
        bfskAvgError = bfskError + bfskAvgError;
    end
    
    ookErrorRate(i) = ookAvgError / testSamples;
    orig_ookErrorRate(i) = orig_ookAvgError / testSamples;
    bpskErrorRate(i) = bpskAvgError / testSamples;
    bfskErrorRate(i) = bfskAvgError / testSamples;

end

% Plot OOK vs DBSK bit error rate
figure(1)
p1 = semilogy(SNR_dB, ookErrorRate,'r-*');
hold on
p2 = semilogy(SNR_dB, bpskErrorRate, 'b-*');
p3 = semilogy(SNR_dB, bfskErrorRate, 'g-*');
p4 = semilogy(SNR_dB, orig_ookErrorRate, 'k-*');

hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1) p3(1) p4(1)],{'Hamming/OOK','Hamming/BPSK',' Hamming/BFSK','Unencoded/OOK'})
xlim([0 50]);
