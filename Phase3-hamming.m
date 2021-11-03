clear all;close all; clc;
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

orig_time = bits/dataRate;      %get the time in seconds
orig_timeScale = 0 : period: orig_time;

%Assume a 6th order Butterworth filter with 0.2 normalised cutoff freq
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
        
%generate data
generatedData = round(randi([0 1], bits,1 ));
generatedData = transpose(generatedData);

%encoding - hamming
encodedHamming = encode(generatedData, codeword_length, message_length ,'hamming/binary');
encodedHamming = transpose(encodedHamming);

%every 160 will change number (either 0/1)
t = dataRate/samplingFreq;
dataSignal = zeros(1, signalLength);
orig_dataSignal = zeros(1, orig_SignalLength);

%for loop cannot start from 0 so need to start from 1. hence, need to -1 to the length
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
ookEnergy = sum(abs(ookSignal).^2);
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
bfskTime = signalLength;
bfskSignalPower = bfskEnergy/bfskTime;
%avgBFSKNoisePower = bfskSignalPower ./ SNR;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;
    
	for j = 1 : testSamples
      
        generatedNoise = randn(1,signalLength);
        orig_generatedNoise = randn(1, orig_SignalLength);
        snrVal = SNR(i);
        % encoded OOK
        receivedOOKSignal = receivedSignal(ookSignal, ookSignalPower, generatedNoise, snrVal);
         
        % unencoded OOK
        orig_receivedOOKSignal = receivedSignal(orig_ookSignal, orig_ookSignalPower, orig_generatedNoise, snrVal);
        
        % BPSK
        receivedBPSKSignal = receivedSignal(bpskSignal, bpskSignalPower, generatedNoise, snrVal);
        
        % BFSK
        receivedBFSKSignal = receivedSignal(bfskSignal, bfskSignalPower, generatedNoise, snrVal);
        
        % Non-coherent Detection method of demodulation
        ookFiltered = ookbpskDemo(receivedOOKSignal,carrierSignal, b, a);
        
        orig_ookDemodulated = orig_receivedOOKSignal .* (2 .* carrierSignal(1:orig_SignalLength));
        orig_ookFiltered = filtfilt(b, a, orig_ookDemodulated);       
        
        bpskFiltered = ookbpskDemo(receivedBPSKSignal,carrierSignal, b, a);
        differenceOfBFSK = bfskDemodulation(receivedBFSKSignal,carrierSignalBFSK1, carrierSignalBFSK2, b, a);
 
        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;
        
        [ookInput, ookOutput] = sampleAndThreshold(ookFiltered, samplingPeriod, avgPower, extended_bits);
        [orig_ookInput, orig_ookOutput] = sampleAndThreshold(orig_ookFiltered, samplingPeriod, avgPower, bits);
        [bpskInput, bpskOutput] = sampleAndThreshold(bpskFiltered, samplingPeriod, 0, extended_bits);
        [bfskInput, bfskOutput] = sampleAndThreshold(differenceOfBFSK, samplingPeriod, 0, extended_bits);

        %Cyclic Decoding
        decoded_OOK = decode(ookOutput,codeword_length, message_length,'hamming/binary');
        decoded_BPSK = decode(bpskOutput,codeword_length, message_length,'hamming/binary');
        decoded_BFSK = decode(bfskOutput,codeword_length, message_length,'hamming/binary');
        
        ookError = biterr(decoded_OOK, generatedData) ./ bits;
        ookAvgError = ookError + ookAvgError;
        
        orig_ookError = biterr(orig_ookOutput, generatedData) ./ bits;
        orig_ookAvgError = orig_ookError + orig_ookAvgError;
        
        bpskError = biterr(decoded_BPSK, generatedData) ./ bits;
        bpskAvgError = bpskError + bpskAvgError;
        
        bfskError = biterr(decoded_BFSK, generatedData) ./ bits;
        bfskAvgError = bfskError + bfskAvgError;
    end
    
	BER_OOK(i) = ookAvgError / testSamples;
    BER_ORIG_OOK(i) = orig_ookAvgError / testSamples;
    BER_BPSK(i) = bpskAvgError / testSamples;
    BER_BFSK(i) = bfskAvgError / testSamples;

end

ookErrorRate = zeros(length(SNR)); %encoded OOK error rate
orig_ookErrorRate = zeros(length(SNR)); %unencoded OOK error rate
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bfskErrorRate = zeros(length(SNR)); %BPSK error rate

% Plot OOK vs DBSK bit error rate
figure(1)
p1 = semilogy(snrDB, BER_OOK,'r-*');
hold on
p2 = semilogy(snrDB, BER_ORIG_OOK, 'k-*');
hold on
p3 = semilogy(snrDB, BER_BPSK, 'b-*');
hold on
p4 = semilogy(snrDB, BER_BFSK, 'g-*');
hold on
legend([p1(1) p2(1) p3(1) p4(1)],{'ENCODED OOK','UNENCODED OOK','BPSK','BFSK'})
xlim([0 50]);
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
hold off
