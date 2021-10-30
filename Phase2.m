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
snrDB = 0:1:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 100;


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
        generatedNoise = randn(1,signalLength);
        
        snrVal = SNR(i);

        receivedOOKSignal = receivedSignal(ookSignal, ookSignalPower, generatedNoise, snrVal);
        receivedBPSKSignal = receivedSignal(bpskSignal, bpskSignalPower, generatedNoise, snrVal);
        receivedBFSKSignal = receivedSignal(bfskSignal, bfskSignalPower, generatedNoise, snrVal);
        
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
       title("Baseband signal oversampled");
        
       subplot(4, 1, 2);
       plot(receivedOOKSignal(1:limit));
       title("OOK modulated signal");

       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 3);
       plot(ookSignal(1:limit));
       title("OOK received with noise");
       
        %Plotting of demodulated signal (mixed and passed through low pass filter
        figure(5);
        subplot(3, 1, 1);
        plot(ookInput(1:limit));
        title("OOK sampled demodulation");

       figure(6);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(ookOutput);
       title("OOK decoded data");
       xlim([0 limit])

       %PLOT BPSK
       figure(7);
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 1);
       plot(dataSignal(1:limit));
       title("Baseband signal oversampled");
        
       subplot(4, 1, 2);
       plot(bpskSignal(1:limit));
       title("BPSK modulated signal");

       subplot(4,1,3);
       plot(receivedBPSKSignal(1:limit));
       title("BPSK received with noise");

       figure(8);
       subplot(3, 1, 1);
       plot(bpskInput(1:limit));
       title("BPSK sampled demodulation");

       figure(9);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bpskOutput);
       title("BPSK decoded data");
       xlim([0 limit])

       %BFSK
       figure(10);
       subplot(4,1,1);
       plot(dataSignal(1:limit));
       title("Baseband signal oversampled");
        

       subplot(4, 1, 2);
       plot(bfskSignal(1:limit));
       title("BFSK modulated signal");

       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 3);
       plot(receivedBFSKSignal(1:limit));
       title("BFSK received with noise");
        

       figure(11);
       subplot(3, 1, 1);
       plot(bfskInput(1:limit));
       title("BFSK sampled demodulation");

       figure(12);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 limit])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bfskOutput);
       title("BFSK decoded data");
       xlim([0 limit])
       
    end
    BER_OOK(i) =  ookAvgError / testSamples;
    BER_BPSK(i) = bpskAvgError / testSamples;
    BER_BFSK(i) = bfskAvgError / testSamples;      

end    

 % Calculate OOK coherent
 coherentOOK = coherent(0, ookNoisePower, amplitude, bits);
 
 % Calculate BPSK coherent
 coherentBPSK = coherent(2, bpskNoisePower, amplitude, bits);

%Calculate BFSK coherent
coherentBFSK = coherent(1, bfskNoisePower, amplitude, bits);


figure('Name','Experimented Data');
plot1 = semilogy(snrDB, BER_OOK, 'b-*');
title('Experimented BER Data');
hold on
plot2 = semilogy (snrDB,BER_BPSK,'r-*');
hold on
plot3 = semilogy (snrDB,BER_BFSK,'g-*');
hold on
legend('OOK','BPSK', 'BFSK');
axis([0 15 10^(-10) 1]);
xlabel('snrDB');
ylabel('BER');
hold off


figure('Name','Theoretical Data');
semilogy (snrDB,coherentOOK,'b-*');
title('Theoretical BER Data');
hold on
semilogy (snrDB,coherentBPSK,'r-*');
hold on
semilogy (snrDB,coherentBFSK,'g-*');
hold on
legend('coherent OOK', 'coherent BPSK', 'coherent BFSK');
axis([0 50 10^(-50) 1]);
xlabel('snrDB');
ylabel('BER');
hold off


        


 








