%Setup
clear all;
close all;
clc;

%carrier frequency:
carrierFreq = 10000; %10kHz for carrier frequency
carrierFreqFSK1 = 10000;        %FSK modulation use - 10kHz
carrierFreqFSK2 = carrierFreqFSK1 * 4;      %FSK modulation use - 40kHz

%Codeword length (n) & Message length (k)
codeword_length = 7;
message_length = 4;
%Special Note: 7/4 seems to be the minimum required

%carrier signal 16 times oversampled:
samplingFreq = 16 * carrierFreq;

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;     %defined from Phase 1
extended_bits = bits*codeword_length/message_length;

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%amplitude for gain
amplitude = 8; 

%% time scale in seconds for ....
time = extended_bits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%Assume a 6th order bandpass filter with cut-off frequency 0.2:
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%carrier signal
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);
carrierSignal_FSK1 = amplitude .* cos(2*pi*carrierFreqFSK1*timeScale);
carrierSignal_FSK2 = amplitude .* cos(2*pi*carrierFreqFSK2*timeScale);

%signal length is for everywhere
signalLength = samplingFreq * extended_bits/dataRate + 1;     %encoded signal length
orig_SignalLength = samplingFreq * bits/dataRate + 1;      %unencoded signal length

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 20;

%generator and syndrome table for cyclic code
%generator polynomial and parity check matrix for cyclic encoding
pol = cyclpoly(codeword_length, message_length);
h = cyclgen(codeword_length, pol);

%syndrome table for cyclic decoding
syndrometable = syndtable(h);

ookErrorRate = zeros(length(SNR)); %encoded OOK error rate
orig_ookErrorRate = zeros(length(SNR)); %unencoded OOK error rate
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bfskErrorRate = zeros(length(SNR)); %BPSK error rate

%generate data
generatedData = round(randi([0 1], bits, 1));
generatedData = transpose(generatedData);

%encoding - cyclic
encodedData = encode(generatedData, codeword_length, message_length, 'cyclic/binary', pol);
encodedData = transpose(encodedData);

%sampling
dataSignal = zeros(1, signalLength);
orig_dataSignal = zeros(1, orig_SignalLength);

for n = 1: signalLength - 1             %encoded signal
    dataSignal(n) = encodedData(ceil(n*dataRate/samplingFreq));
end

for n = 1 : orig_SignalLength - 1          %original signal
    orig_dataSignal(n) = generatedData(ceil(n*dataRate/samplingFreq));
end

%% Why overwrite?
dataSignal(signalLength) = dataSignal(signalLength - 1);
orig_dataSignal(orig_SignalLength) = orig_dataSignal(orig_SignalLength - 1);

%==== OOK ====%
%encoded signal
ookSignal = carrierSignal .* dataSignal;
ookSignalPower = (norm(ookSignal)^2)/signalLength;
%ookNoisePower = ookSignalPower ./ SNR;

%unencoded signal
orig_ookSignal = carrierSignal(1:orig_SignalLength) .* orig_dataSignal;
orig_ookSignalPower = (norm(orig_ookSignal)^2)/orig_SignalLength;
%orig_ookNoisePower = orig_ookSignalPower ./ SNR;

%==== BPSK ====%
bpskSourceSignal = dataSignal .* 2 - 1; %*2 bc of negative 1 and 1 bc ook is 0 to 1 so don't need *2
bpskSignal = carrierSignal .* bpskSourceSignal;
bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
%bpskNoisePower = bpskSignalPower ./ SNR;

%==== BFSK ====%
bfskSourceSignal_high = carrierSignal_FSK1 .* (dataSignal == 1);
bfskSourceSignal_low = carrierSignal_FSK2 .* (dataSignal == 0);
bfskSignal = bfskSourceSignal_low + bfskSourceSignal_high;
bfskSignalPower = (norm(bfskSignal)^2)/signalLength;
%bfskNoisePower = bfskSignalPower ./ SNR;

% For each value of SNR, test of 100 samples
for i = 1: length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;
    result = zeros(1, testSamples);

    for j = 1 : testSamples
        %encoded OOK
        avgOOKNoisePower = ookSignalPower ./ SNR(i);
        ookNoise = sqrt(avgOOKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedOOKSignal = ookSignal + ookNoise;

        %unencoded OOK
        orig_avgOOKNoisePower = orig_ookSignalPower ./ SNR(i);
        orig_ookNoise = sqrt(orig_avgOOKNoisePower) .* transpose(randi([0 1], orig_SignalLength, 1));
        orig_receivedOOKSignal = orig_ookSignal + orig_ookNoise;

        %BPSK
        avgBPSKNoisePower = bpskSignalPower ./ SNR(i);
        bpskNoise = sqrt(avgBPSKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedBPSKSignal = bpskSignal + bpskNoise;

        %BFSK
        avgBFSKNoisePower = bfskSignalPower ./ SNR(i);
        bfskNoise = sqrt(avgBFSKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedBFSKSignal = bfskSignal + bfskNoise;

        %Non-coherent Detection method of demodulation
        ookDemodulated = receivedOOKSignal .* (2 .* carrierSignal);
        ookFiltered = filtfilt(b, a, ookDemodulated);

        orig_ookDemodulated = orig_receivedOOKSignal .* (2 .* carrierSignal(1:orig_SignalLength));
        orig_ookFiltered = filtfilt(b, a, orig_ookDemodulated);

        bpskDemodulated = receivedBPSKSignal .* (2 .* carrierSignal);
        bpskFiltered = filtfilt(b, a, bpskDemodulated);

        bfskDemodulated_high = receivedBFSKSignal .* (2 .* carrierSignal_FSK1);
        bfskFiltered_high = filtfilt(b, a, bfskDemodulated_high);
        bfskDemodulated_low = receivedBFSKSignal .* (2 .* carrierSignal_FSK2);
        bfskFiltered_low = filtfilt(b, a, bfskDemodulated_low);
        bfskFiltered = bfskFiltered_high - bfskFiltered_low;

        %sampling period for demodulation
        %% Question: isn't samplingPeriod = samplingRate (line 26)?
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;

        [ookInput, ookOutput] = samplingforthreshold(ookFiltered, samplingPeriod, avgPower, extended_bits);
        [orig_ookInput, orig_ookOutput] = samplingforthreshold(orig_ookFiltered, samplingPeriod, avgPower, bits);    %note the change to unencoded bits
        [bpskInput, bpskOutput] = samplingforthreshold(bpskFiltered, samplingPeriod, 0, extended_bits); %threshold is 0 
        [bfskInput, bfskOutput] = samplingforthreshold(bfskFiltered, samplingPeriod, 0, extended_bits);

        %Cyclic Decoding
        decoded_OOK = decode(ookOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        decoded_BPSK = decode(bpskOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        decoded_BFSK = decode(bfskOutput, codeword_length, message_length, 'cyclic/binary', pol, syndrometable);
        
        ookError = biterr(decoded_OOK, generatedData) ./ bits;
        orig_ookError = biterr(orig_ookOutput, generatedData) ./ bits;
        bpskError = biterr(decoded_BPSK, generatedData) ./ bits;
        bfskError = biterr(decoded_BFSK, generatedData) ./ bits;

        ookAvgError = ookAvgError + ookError;
        orig_ookAvgError = orig_ookAvgError + orig_ookError;
        bpskAvgError = bpskAvgError + bpskError;
        bfskAvgError = bfskAvgError + bfskError;
    end

    ookErrorRate(i) = ookAvgError / testSamples;
    orig_ookErrorRate(i) = orig_ookAvgError / testSamples;
    bpskErrorRate(i) = bpskAvgError / testSamples;
    bfskErrorRate(i) = bfskAvgError / testSamples;
    
    %Plot for SNR @ 5dB?
    %if(i == 5)
        %% To continue

    %end
end


%OOK vs BPSK vs BFSK bit error rate
figure(1)
p1 = semilogy(snrDB, ookErrorRate, 'r-*');
hold on
p2 = semilogy(snrDB, orig_ookErrorRate, 'k-*');
p3 = semilogy(snrDB, bpskErrorRate, 'b-*');
p4 = semilogy(snrDB, bfskErrorRate, 'g-*');
hold off

ylabel ('BER - Bit Error Rate');
xlabel ('SNR(dB) - Signal to Noise Ratio in dB');
legend([p1(1) p2(1) p3(1) p4(1)], {'Cyclic - Encoded OOK', 'Cyclic - Unencoded OOK', 'Cyclic - BPSK', 'Cyclic - BFSK'})
%legend([p1(1) p2(1)], {'Cyclic - Encoded OOK', 'Cyclic - Unencoded OOK'})
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
