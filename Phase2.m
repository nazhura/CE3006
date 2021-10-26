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

%number of test per samples
testSamples = 100;

ookErrorRate = zeros(length(SNR)); %OOK error rate
bpskErrorRate = zeros(length(SNR)); %BPSK error rate
bfskErrorRate = zeros(length(SNR)); %BFSK error rate

%generate data
generatedData = round(randi([0 1], bits, 1));
generatedData = transpose(generatedData);

dataSignal = zeros(1, signalLength);

for n = 1: signalLength - 1
    dataSignal(n) = generatedData(ceil(n*dataRate/samplingFreq));
end

dataSignal(signalLength) = dataSignal(signalLength - 1);

%==== OOK ====%
ookSignal = carrierSignal .* dataSignal;
ookSignalPower = (norm(ookSignal)^2)/signalLength;
ookNoisePower = ookSignalPower ./ SNR;

%==== BPSK ====%
bpskSourceSignal = dataSignal .* 2 - 1; %*2 bc of negative 1 and 1 bc ook is 0 to 1 so don't need *2 (project manual got say)
bpskSignal = carrierSignal .* bpskSourceSignal; % this too
bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
bpskNoisePower = bpskSignalPower ./ SNR;

%==== BFSK ===%
carrierFreqBFSK1 = 50000; %first frequency is 50K
carrierFreqBFSK2 = 10000; %second frequency is 10K

%generate carrier signal for BFSK:
carrierSignalBFSK1 = amplitude .* cos(2*pi*carrierFreqBFSK1*timeScale);
carrierSignalBFSK2 = amplitude .* cos(2*pi*carrierFreqBFSK2*timeScale);

%BFSK modulation:

%bfskSignal1 = carrierSignalBFSK1 * firstSignal; %signal = 1
%bfskSignal2 = carrierSignalBFSK2 * secondSignal; %signal = 0
%^ for our understanding delete when submitting

bfskSignal = ((carrierSignalBFSK1 .* (dataSignal == 1)) + (carrierSignalBFSK2 .* (dataSignal == 0)));
bfskSignalPower = (norm(bfskSignal)^2)/signalLength;
bfskNoisePower = bfskSignalPower ./ SNR;


for i = 1 : length(SNR)
    ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;

    for j = 1 : testSamples
        avgOOKNoisePower = ookSignalPower ./ SNR(i);
        ookNoise = sqrt(avgOOKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedOOKSignal = ookSignal + ookNoise;

        avgBPSKNoisePower = bpskSignalPower ./ SNR(i);
        bpskNoise = sqrt(avgBPSKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedBPSKSignal = bpskSignal + bpskNoise;

        avgBFSKNoisePower = bfskSignalPower ./ SNR(i);
        bfskNoise = sqrt(avgBFSKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receviedBFSKSignal = bfskSignal + bfskNoise;


        %demodulation
        ookDemodulated = receivedOOKSignal .* (2 .* carrierSignal);
        ookFiltered = filtfilt(b, a, ookDemodulated);

        bpskDemodulated = receivedBPSKSignal .* (2 .* carrierSignal);
        bpskFiltered = filtfilt(b, a, bpskDemodulated);

        bfskDemodulated1 = receviedBFSKSignal .* (2 .* carrierSignalBFSK1);
        bfskDemodulated2 = receviedBFSKSignal .* (2 .* carrierSignalBFSK2);
        bfskFiltered1 = filtfilt(b, a, bfskDemodulated1);
        bfskFiltered2 = filtfilt(b, a, bfskDemodulated2);
        differenceOfBFSK = bfskFiltered1 - bfskFiltered2;

        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;
        
        [ookInput, ookOutput] = samplingforthreshold(ookFiltered, samplingPeriod, avgPower, bits);
        [bpskInput, bpskOutput] = samplingforthreshold(bpskFiltered, samplingPeriod, 0, bits); %threshold is 0 
        [bfskInput, bfskOutput] = samplingforthreshold(differenceOfBFSK, samplingPeriod, 0, bits); %threshold is 0


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
       figure(2)
       plot(dataSignal, 'b');
       title("Original Data")
       xlim([0 1024])

       %PLOT OOK
       figure(3);
       subplot(4,1,1);
       plot(dataSignal(1:1024));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receivedOOKSignal(1:1024));
       title("OOK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(ookSignal(1:1024));
       title("OOK modulated signal (Snippet)");

       
        %Plotting of demodulated signal (mixed and passed through low pass filter
        figure(5);
        subplot(3, 1, 1);
        plot(ookInput(1:1024));
        title("OOK demodulated and sampled");

       figure(6);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 1024])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(ookOutput);
       title("OOK Decoded Data");
       xlim([0 1024])

       %PLOT BPSK
       figure(7);
       subplot(4,1,1);
       plot(dataSignal(1:1024));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receivedBPSKSignal(1:1024));
       title("BPSK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(bpskSignal(1:1024));
       title("BPSK modulated signal (Snippet)");

       figure(8);
       subplot(3, 1, 1);
       plot(bpskInput(1:1024));
       title("BPSK demodulated and sampled");

       figure(9);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 1024])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bpskOutput);
       title("BPSK Decoded Data");
       xlim([0 1024])

       %BFSK
       figure(10);
       subplot(4,1,1);
       plot(dataSignal(1:1024));
       title("Baseband oversampled signal (snippet)");
        
       %Plotting of Received signal (corrupted with noise)
       subplot(4, 1, 2);
       plot(receviedBFSKSignal(1:1024));
       title("BFSK received with noise (Snippet)");
        
       subplot(4, 1, 3);
       plot(bfskSignal(1:1024));
       title("BFSK modulated signal (Snippet)");

       figure(11);
       subplot(3, 1, 1);
       plot(bfskInput(1:1024));
       title("BFSK demodulated and sampled");

       figure(12);
       subplot(4, 1, 1);
       plot(generatedData);
       title("Original Data")
       xlim([0 1024])
       ylim([0 1])

       subplot(4, 1, 2);
       plot(bfskOutput);
       title("BFSK Decoded Data");
       xlim([0 1024])
       

    end

end    

% Calculate OOK coherent
e1OOK = (1 / 2) * amplitude^2 / bits;
e0OOk = 0;
ebOOK = (1 / 2) * (e1OOK + e0OOk);
noOOK = ookNoisePower ./ bits ./ 2;
coherentOOK = (1 / 2) .* erfc(sqrt(ebOOK ./ (2 .* noOOK)));

% Calculate BPSK coherent
e1BPSK = (1 / 2) * amplitude^2 / bits;
e0BPSK = (1 / 2) * amplitude^2 / bits;
ebBPSK = (1 / 2) * (e1BPSK + e0BPSK);
noBPSK = bpskNoisePower ./ bits ./ 2;
coherentBPSK = (1 / 2) .* erfc(sqrt(ebBPSK ./ (2 .* noBPSK)));

% Calculate BFSK coherent
e1BFSK = (1 / 2) * amplitude^2 / bits;
e0BFSK = (1 / 2) * amplitude^2 / bits;
ebBFSK = (1 / 2) * (e1BFSK + e0BFSK);
noBFSK = bfskNoisePower ./ bits ./ 2;
coherentBFSK = (1 / 2) .* erfc(sqrt(ebBFSK ./ (2 .* noBFSK)));
        

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
 


%need plot semilogy for all:
