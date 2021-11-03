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

%Self-defined: Codeword length (n) & Message length (k)
codeword_length = 7;
message_length = 4;

%carrier signal 16 times oversampled:
samplingFreq = 16* carrierFreq; %sampling frequency is 16 times the carrier frequency

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;
extended_bits = bits*codeword_length/message_length;

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%amplitude for
amplitude = 8;

%timescale in seconds for ....
time = extended_bits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%assume a 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%carrier signal
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);

%signal length is for everywhere
signalLength = samplingFreq*extended_bits/dataRate + 1;
orig_SignalLength = samplingFreq* bits/dataRate + 1;

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 10;

ookErrorRate = zeros(length(SNR));
ookErrorRate_org = zeros(length(SNR));
bpskErrorRate = zeros(length(SNR));
bfskErrorRate = zeros(length(SNR));
        
%generate noise - Hamming
generatedData = round(randi([0 1], bits,1 ));
generatedData = transpose(generatedData);

%encoding - hamming
encodedHamming = encode(generatedData, codeword_length, message_length ,'hamming/binary');
encodedHamming = transpose(encodedHamming);

%sampling
t = dataRate/samplingFreq;
dataSignal = zeros(1, signalLength);
orig_dataSignal = zeros(1, orig_SignalLength);
for n = 1: signalLength - 1
    dataSignal(n) = encodedHamming(ceil(n*t));
end
for n = 1: orig_SignalLength - 1
    orig_dataSignal(n) = generatedData(ceil(n*t));
end
dataSignal(signalLength) = dataSignal(signalLength - 1);
orig_dataSignal(orig_SignalLength) = orig_dataSignal(orig_SignalLength - 1);

%==== OOK ====%
%encoded signal
ookSignal = carrierSignal .* dataSignal;
ookEnergy = sum(abs(norm(ookSignal).^2);
ookTime = signalLength;
ookSignalPower = ookEnergy/ookTime; %Power = Energy/Time

%unencoded signal
orig_ookSignal = carrierSignal(1:orig_SignalLength) .* orig_dataSignal;
orig_ookEnergy = sum(abs(orig_ookSignal).^2);
orig_ookTime = orig_SignalLength;
orig_ookSignalPower = orig_ookEnergy/orig_ookTime;

%==== BPSK ====%
bpskSourceSignal = 2.*(dataSignal-0.5);
bpskSignal = carrierSignal .* bpskSourceSignal;
bpskEnergy = sum(abs(bpskSignal).^2);
bpskTime = signalLength;
bpskSignalPower = bpskEnergy/bpskTime;

%==== BFSK ====%
carrierFreqBFSK1 = 50000;
carrierFreqBFSK2 = 10000;

carrierSignalBFSK1 = amplitude .* cos(2*pi*carrierFreqBFSK1*timeScale);
carrierSignalBFSK2 = amplitude .* cos(2*pi*carrierFreqBFSK2*timeScale);

%BFSK modulation:
dataSignal2 = mod(dataSignal + 1, 2); %signal = 1

bfskSignal = carrierSignalBFSK1 .* dataSignal + carrierSignalBFSK2 .* dataSignal2;
bfskEnergy = sum(abs(bfskSignal).^2);
bfsktime = signalLength;
bfskSignalPower = bfskEnergy/bfskTime;
%avgBFSKNoisePower = bfskSignalPower ./ SNR;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;
    
	for j = 1 : testSamples
      
        % encoded OOK
        avgOOKNoisePower = ookSignalPower ./SNR(i);
        ookNoise = sqrt(avgOOKNoisePower) .*randn(1, signalLength);
        receivedOOKSignal = ookSignal + ookNoise;
        
        % unencoded OOK
        orig_avgOOKNoisePower = orig_ookSignalPower ./ SNR(i);
        orig_ookNoise =  sqrt(orig_avgOOKNoisePower) .* randn(1, orig_SignalLength);
        orig_receivedOOKSignal = orig_ookSignal + orig_ookNoise;
        
        % BPSK
        avgBPSKNoisePower = bpskSignalPower ./SNR(i);
        bpskNoise = sqrt(avgBPSKNoisePower) .* randn(1, signalLength);
        receivedBPSKSignal = bpskSignal + bpskNoise;
        
        % BFSK
        avgBFSKNoisePower = bfskSignalPower ./SNR(i);
        bfskNoise = sqrt(avgBFSKNoisePower) .*randn(1, signalLength);
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
        
        [ookInput, ookOutput] = samplingforthreshold(ookFiltered, samplingPeriod, avgPower, extended_bits);
        [orig_ookInput, orig_ookOutput] = samplingforthreshold(orig_ookFiltered, samplingPeriod, avgPower, bits);
        [bpskInput, bpskOutput] = samplingforthreshold(bpskFiltered, samplingPeriod, 0, extended_bits);
        [bfskInput, bfskOutput] = samplingforthreshold(bfskFiltered, samplingPeriod, 0, extended_bits);

        %Cyclic Decoding
        decoded_OOK = decode(ookOutput,codeword_length, message_length,'hamming/binary');
        decoded_BPSK = decode(bpskOutput,codeword_length, message_length,'hamming/binary');
        decoded_BFSK = decode(bfskOutput,codeword_length, message_length,'hamming/binary');
        
        ookError = biterr(decoded_OOK, generatedData) ./ bits;
        orig_ookError = biterr(orig_ookOutput, generatedData) ./ bits;
        bpskError = biterr(decoded_BPSK, generatedData) ./ bits;
        bfskError = biterr(decoded_BFSK, generatedData) ./ bits;
        
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
p1 = semilogy(snrDB, ookErrorRate,'r-*');
hold on
p2 = semilogy(snrDB, bpskErrorRate, 'b-*');
p3 = semilogy(snrDB, bfskErrorRate, 'g-*');
p4 = semilogy(snrDB, orig_ookErrorRate, 'k-*');

hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1) p3(1) p4(1)],{'Hamming/OOK','Hamming/BPSK',' Hamming/BFSK','Unencoded/OOK'})
xlim([0 50]);

function [input, output] = samplingforthreshold(filter, period, threshold, bits)
    input = zeros(1, bits);
    output = input;
    avgTime = period / 2;

    for i = 1: bits
        input(i) = filter((2 * i - 1) * avgTime);
        if(input(i) > threshold)
            output(i) = 1;
        else
            output(i) = 0;
        end
    end
end