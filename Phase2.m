clear all; close all; clc;

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

%define amplitude
amplitude = 8; 

%time scale in seconds for ....
time = bits/dataRate; %get the time in seconds
period = 1/samplingFreq; %period = 1/frequency
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

%number of test per samples
testSamples = 10;


%generate data
generatedData = randi([0 1], bits, 1);
generatedData = transpose(generatedData);


%every 160 will change number (either 0/1)
t = dataRate / samplingFreq;
dataSignal = zeros(1, signalLength);

for n = 1: signalLength - 1 %for loop cannot start from 0 so need to start from 1. hence, need to -1 to the length
    dataSignal(n) = generatedData(ceil(n*t));
end

dataSignal(signalLength) = dataSignal(signalLength - 1);
%taking last dataSignal (163841) to be the same as 163840 

%==== OOK ====%
ookSignal = carrierSignal .* dataSignal;
ookEnergy = sum(abs(ookSignal).^2);
ookTime = signalLength;

ookSignalPower = ookEnergy/ookTime; %Power = Energy/Time
ookNoisePower = ookSignalPower ./ SNR;





%==== BPSK ====%
% Convert Input Binary Data to +1 and -1
% If Input Binary Data = 1, 2 * (1 - 0.5) = 1
% If input Binary Data = 0, 2 * (0 - 0.5) = -1
bpskSourceSignal = 2 .* (dataSignal-0.5); 

bpskSignal = carrierSignal .* bpskSourceSignal; 
bpskEnergy = sum(abs(bpskSignal).^2);
bpskTime = signalLength;

bpskSignalPower = bpskEnergy/bpskTime;
bpskNoisePower = bpskSignalPower ./ SNR;


%==== BFSK ===%
carrierFreqBFSK1 = 50000; %first frequency is 50K
carrierFreqBFSK2 = 10000; %second frequency is 10K

%generate carrier signal for BFSK:
carrierSignalBFSK1 = amplitude .* cos(2*pi*carrierFreqBFSK1*timeScale);
carrierSignalBFSK2 = amplitude .* cos(2*pi*carrierFreqBFSK2*timeScale);

%BFSK modulation:
dataSignal2 = mod(dataSignal + 1, 2); %signal = 1
%Getting 0 for carrierSignal2:
%(1 + 1 = 2) % 2 = 0
%(0 + 1 = 1) % 2 = 1

bfskSignal = carrierSignalBFSK1 .* dataSignal + carrierSignalBFSK2 .* dataSignal2;
bfskEnergy = sum(abs(bfskSignal).^2);
bfskTime = signalLength;

bfskSignalPower = bfskEnergy/bfskTime;
bfskNoisePower = bfskSignalPower ./ SNR;


for i = 1 : length(SNR)
    ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;

    for j = 1 : testSamples

        %generate data
        generatedSignal = randi([0 1], signalLength, 1);
        generatedSignal = transpose(generatedSignal);

        snrVal = SNR(i);

        receivedOOKSignal = receivedSignal(ookSignal, ookSignalPower, generatedSignal, snrVal);
        receivedBPSKSignal = receivedSignal(bpskSignal, bpskSignalPower, generatedSignal, snrVal);
        receivedBFSKSignal = receivedSignal(bfskSignal, bfskSignalPower, generatedSignal, snrVal);
        
        %demodulation

        ookFiltered = ookbpskDemo(receivedOOKSignal,carrierSignal, b, a);
        bpskFiltered = ookbpskDemo(receivedBPSKSignal,carrierSignal, b, a);
        differenceOfBFSK = bfskDemodulation(receivedBFSKSignal,carrierSignalBFSK1, carrierSignalBFSK2, b, a);

        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;

        [ookInput, ookOutput] = sampleAndThreshold(ookFiltered, samplingPeriod, avgPower, bits);
        [bpskInput, bpskOutput] = sampleAndThreshold(bpskFiltered, samplingPeriod, 0, bits); %threshold is 0 
        [bfskInput, bfskOutput] = sampleAndThreshold(differenceOfBFSK, samplingPeriod, 0, bits); %threshold is 0

        ookError = 0;
        bpskError = 0;
        bfskError = 0;

        for k = 1: bits
            if(ookOutput(k) ~= generatedData(k))
                ookError = ookError + 1;
            end
             if(bpskOutput(k) ~= generatedData(k))
                bpskError = bpskError + 1;
             end
             if(bfskOutput(k) ~= generatedData(k))
                 bfskError = bfskError + 1;
             end

        end

        ookError = ookError./bits;
        ookAvgError = ookError + ookAvgError;

        bpskError = bpskError./bits;
        bpskAvgError = bpskError + bpskAvgError;

        bfskError = bfskError./bits;
        bfskAvgError = bfskError + bfskAvgError;
        
    end

    %Plot the 5db SNR signals
    if (snrDB(i) == 5)
    %Plot of original data with respect to time
       limit = 1024;

       figure(1)
       plot(generatedData, 'b');
       title("Generated Data")
       xlim([0 limit])

       figure(2)
       plot(dataSignal, 'b');
       title("Original Data")
       xlim([0 limit])

       %PLOT OOK
       figure(3);
       subplot(4,1,1);
       plot(dataSignal(1:limit));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receivedOOKSignal(1:limit));
       title("OOK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(ookSignal(1:limit));
       title("OOK modulated signal (Snippet)");

       
        %Plotting of demodulated signal (mixed and passed through low pass filter
        figure(5);
        subplot(3, 1, 1);
        plot(ookInput(1:limit));
        title("OOK demodulated and sampled");

       figure(6);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(ookOutput);
       title("OOK Decoded Data");
       xlim([0 limit])

       %PLOT BPSK
       figure(7);
       subplot(4,1,1);
       plot(dataSignal(1:limit));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receivedBPSKSignal(1:limit));
       title("BPSK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(bpskSignal(1:limit));
       title("BPSK modulated signal (Snippet)");

       figure(8);
       subplot(3, 1, 1);
       plot(bpskInput(1:limit));
       title("BPSK demodulated and sampled");

       figure(9);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bpskOutput);
       title("BPSK Decoded Data");
       xlim([0 limit])

       %BFSK
       figure(10);
       subplot(4,1,1);
       plot(dataSignal(1:limit));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receivedBFSKSignal(1:limit));
       title("BFSK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(bfskSignal(1:limit));
       title("BFSK modulated signal (Snippet)");

       figure(11);
       subplot(3, 1, 1);
       plot(bfskInput(1:limit));
       title("BFSK demodulated and sampled");

       figure(12);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bfskOutput);
       title("BFSK Decoded Data");
       xlim([0 limit])
       
    end
    BER_OOK(i) =  ookAvgError / testSamples;
    BER_BPSK(i) = bpskAvgError / testSamples;
    BER_BFSK(i) = bfskAvgError / testSamples;      

end    

% Calculate OOK coherent
coherentOOK = coherent(0, ookNoisePower, amplitude, bits);

% Calculate BPSK coherent
coherentBPSK = coherent(1, bpskNoisePower, amplitude, bits);

% Calculate BFSK coherent
coherentBFSK = coherent(1, bfskNoisePower, amplitude, bits);

%need plot semilogy for all:
ookErrorRate = zeros(length(SNR)); %OOK error rate
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bfskErrorRate = zeros(length(SNR)); %BFSK error rate

figure(13);
figure('Name','Measured Data');
semilogy(snrDB, BER_OOK, 'b-*');
title('BER against SNR');
hold on
semilogy (snrDB,BER_BPSK,'r-*');
hold on
semilogy (snrDB,BER_BFSK,'g-*');
hold on
semilogy (snrDB,coherentOOK,'m-*');
hold on
semilogy (snrDB,coherentBPSK,'k-*');
hold on
semilogy (snrDB,coherentBFSK,'c-*');
hold on
legend('OOK','BPSK', 'BFSK');
axis([0 30 10^(-5) 1]);
xlabel('snrDB');
ylabel('BER');
hold off



        


 








