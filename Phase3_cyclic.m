clear all; close all; clc;

%carrier frequency:
carrierFreq = 10000; %10kHz for carrier frequency

%Self-defined: Codeword length (n) & Message length (k)
codeword_length = 6;
message_length = 4;
%Some possible combinations: 3/1, 7/4, 15/11, 31/26

%carrier signal 16 times oversampled:
samplingFreq = 16 * carrierFreq;

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;     %defined from Phase 1
extended_bits = bits*codeword_length/message_length;    %for unencoded bit length

samplingRate = samplingFreq / dataRate;

%amplitude for gain
amplitude = 8; 

% time scale
time = extended_bits/dataRate;      %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

orig_time = bits/dataRate;      %get the time in seconds
orig_timeScale = 0 : period: orig_time;

%Assume a 6th order bandpass filter with cut-off frequency 0.2:
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%carrier signal
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);
orig_carrierSignal = amplitude .* cos(2*pi*carrierFreq*orig_timeScale);

%signal length is for everywhere
signalLength = samplingFreq * extended_bits/dataRate + 1;     %encoded signal length
orig_SignalLength = samplingFreq * bits/dataRate + 1;      %unencoded signal length

%SNR
snrDB = 0:0.5:15;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 100;

%For Cyclic code: g (poly), h (h) and s (syndrometable)
%generator polynomial and parity check matrix for cyclic encoding
pol = cyclpoly(codeword_length, message_length);
h = cyclgen(codeword_length, pol);

%syndrome table for cyclic decoding
syndrometable = syndtable(h);

%generate data
generatedData = randi([0 1], bits, 1);
generatedData = transpose(generatedData);

%encoding - cyclic
encodedData = encode(generatedData, codeword_length, message_length, 'cyclic/binary', pol);
encodedData = transpose(encodedData);

%every 160 will change number (either 0/1)
t = dataRate / samplingFreq;
dataSignal = zeros(1, signalLength);
orig_dataSignal = zeros(1, orig_SignalLength);

%for loop cannot start from 0 so need to start from 1. hence, need to -1 to the length
for n = 1: signalLength - 1             %encoded signal
    dataSignal(n) = encodedData(ceil(n*t));
end

for n = 1 : orig_SignalLength - 1          %original signal
    orig_dataSignal(n) = generatedData(ceil(n*t));
end

dataSignal(signalLength) = dataSignal(signalLength - 1);
orig_dataSignal(orig_SignalLength) = orig_dataSignal(orig_SignalLength - 1);

%==== OOK ====%
%encoded signal
ookSignal = carrierSignal .* dataSignal;
ookEnergy = sum(abs(ookSignal).^2);
ookTime = signalLength;
ookSignalPower = ookEnergy/ookTime; %Power = Energy/Time

%unencoded signal
orig_ookSignal = orig_carrierSignal .* orig_dataSignal;
orig_ookEnergy = sum(abs(orig_ookSignal).^2);
orig_ookTime = orig_SignalLength;
orig_ookSignalPower = orig_ookEnergy/orig_ookTime;

%==== BPSK ====%
% Convert Input Binary Data to +1 and -1
% If Input Binary Data = 1, 2 * (1 - 0.5) = 1
% If input Binary Data = 0, 2 * (0 - 0.5) = -1
bpskSourceSignal = 2 .* (dataSignal-0.5);

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
dataSignal2 = mod(dataSignal + 1, 2);

bfskSignal = carrierSignalBFSK1 .* dataSignal + carrierSignalBFSK2 .* dataSignal2;
bfskEnergy = sum(abs(bfskSignal).^2);
bfskTime = signalLength;

bfskSignalPower = bfskEnergy/bfskTime;

% For each value of SNR, test of 100 samples
for i = 1: length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;

    for j = 1 : testSamples
        
        generatedNoise = randn(1,signalLength);
        orig_generatedNoise = randn(1, orig_SignalLength);
        snrVal = SNR(i);

        receivedOOKSignal = receivedSignal(ookSignal, ookSignalPower, generatedNoise, snrVal); %encoded OOK
        orig_receivedOOKSignal = receivedSignal(orig_ookSignal, orig_ookSignalPower, orig_generatedNoise, snrVal); %unencoded OOK
        receivedBPSKSignal = receivedSignal(bpskSignal, bpskSignalPower, generatedNoise, snrVal); %encoded BPSK
        receivedBFSKSignal = receivedSignal(bfskSignal, bfskSignalPower, generatedNoise, snrVal); %encoded BFSK

        %demodulation
        ookFiltered = ookbpskDemo(receivedOOKSignal,carrierSignal, b, a);
        orig_ookFiltered = ookbpskDemo(orig_receivedOOKSignal,orig_carrierSignal, b, a); %note the change for unencoded bits
        bpskFiltered = ookbpskDemo(receivedBPSKSignal,carrierSignal, b, a);
        differenceOfBFSK = bfskDemodulation(receivedBFSKSignal,carrierSignalBFSK1, carrierSignalBFSK2, b, a);

        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;

        [ookInput, ookOutput] = sampleAndThreshold(ookFiltered, samplingPeriod, avgPower, extended_bits);
        [orig_ookInput, orig_ookOutput] = sampleAndThreshold(orig_ookFiltered, samplingPeriod, avgPower, bits);    %note the change for unencoded bits
        [bpskInput, bpskOutput] = sampleAndThreshold(bpskFiltered, samplingPeriod, 0, extended_bits); %threshold is 0, not Aˆ2/2 
        [bfskInput, bfskOutput] = sampleAndThreshold(differenceOfBFSK, samplingPeriod, 0, extended_bits);

        %Cyclic Decoding
        decoded_OOK = decode(ookOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        decoded_BPSK = decode(bpskOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        decoded_BFSK = decode(bfskOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        
        ookError = 0;                                           
        orig_ookError = 0;
        bpskError = 0;
        bfskError = 0;

        for k = 1: bits
            if(decoded_OOK(k) ~= generatedData(k))
                ookError = ookError + 1;
            end
            if(orig_ookOutput(k) ~= generatedData(k))
                orig_ookError = orig_ookError+1;
            end
            if(decoded_BPSK(k) ~= generatedData(k))
               bpskError = bpskError + 1;
            end
            if(decoded_BFSK(k) ~= generatedData(k))
               bfskError = bfskError + 1;
            end

        end

        ookError = ookError ./bits;
        ookAvgError = ookAvgError + ookError;
        orig_ookError = orig_ookError ./bits;
        orig_ookAvgError = orig_ookAvgError + orig_ookError;
        bpskError = bpskError ./bits;
        bpskAvgError = bpskAvgError + bpskError;
        bfskError = bfskError ./bits;
        bfskAvgError = bfskAvgError + bfskError;
    end
    %Removed plot
    BER_OOK(i) = ookAvgError / testSamples;
    BER_ORIG_OOK(i) = orig_ookAvgError / testSamples;
    BER_BPSK(i) = bpskAvgError / testSamples;
    BER_BFSK(i) = bfskAvgError / testSamples;
end

%Removed theoretical BER
ookErrorRate = zeros(length(SNR)); %encoded OOK error rate
orig_ookErrorRate = zeros(length(SNR)); %unencoded OOK error rate
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bfskErrorRate = zeros(length(SNR)); %BPSK error rate

figure('Name','Measured Data');
title('BER against SNR');
semilogy(snrDB, BER_OOK, 'b-*');
hold on
semilogy(snrDB, BER_ORIG_OOK, 'k-*');
hold on
semilogy (snrDB,BER_BPSK,'r-*');
hold on
semilogy (snrDB,BER_BFSK,'g-*');
hold on
legend('ENCODED OOK','UNENCODED OOK', 'BPSK', 'BFSK');
axis([0 15 10^(-7) 1]);
xlabel('Signal-to-Noise Ratio (in dB)');
ylabel('Bit Error Rate (BER)');
hold off
