%Setup
clear all;
close all;
clc;

%carrier frequency:
carrierFreq = 10000; %10kHz for carrier frequency
carrierFreqFSK1 = 10000;        %FSK modulation use - 10kHz
carrierFreqFSK2 = carrierFreqFSK1 * 4;      %FSK modulation use - 40kHz

%Self-defined: Codeword length (n) & Message length (k)
codeword_length = 6;
message_length = 4;

%carrier signal 16 times oversampled:
samplingFreq = 16 * carrierFreq; %sampling frequency is 16 times the carrier frequency

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;    %defined from Phase 1
msgbits = bits*codeword_length/message_length;  %actual message bit, excluding parity bits

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%amplitude for 
amplitude = 8; 

%% time scale in seconds for ....
time = msgbits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%Assume a 6th order bandpass filter with cut-off frequency 0.2:
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%carrier signal
carrierSignal = amplitude .* cos(2*pi*carrierFreq*timeScale);
carrierSignal_FSK1 = amplitude .* cos(2*pi*carrierFreqFSK1*timeScale);
carrierSignal_FSK2 = amplitude .* cos(2*pi*carrierFreqFSK2*timeScale);

%signal length is for everywhere
signalLength = samplingFreq * msgbits/dataRate + 1;     %encoded signal length
actualSignalLength = samplingFreq * bits/dataRate + 1;      %unencoded signal length

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 100;

%% What is going on??
%generator and syndrome table for cyclic code
%generator polynomial and parity check matrix for cyclic encoding
polynom = cyclpoly(codeword_length, message_length);
parmat = cyclgen(codeword_length, polynom);

%syndrome table for cyclic decoding
table = syndtable(parmat);

ookErrorRate = zeros(length(SNR)); %OOK error rate
ookErrorRate_org = zeros(length(SNR)); %OOK error rate - unencoded
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bpskErrorRate_org = zeros(length(SNR)); %BPSK error rate

%generate data
generatedData = round(randi([0 1], bits, 1));

%encoding - cyclic
encodedData = encode(generatedData, codeword_length, message_length, 'cyclic/binary', polynom);

%sampling
dataSignal = zeros(1, signalLength);
actualDataSignal = zeros(1, actualSignalLength);

for n = 1: signalLength - 1             %encoded signal
    dataSignal(n) = encodedData(ceil(n*dataRate/samplingFreq));
end

for n = 1 : actualSignalLength - 1          %original signal
    actualDataSignal(n) = generatedData(ceil(n*dataRate/samplingFreq));
end


%% to continue
