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


for i = 1 : length(SNR)
    ookAverageError = 0;

    for j = 1 : testSamples
        avgOOKNoisePower = ookSignalPower ./ SNR(i);
        ookNoise = sqrt(avgOOKNoisePower) .* transpose(randi([0 1], signalLength, 1));
        receivedOOKSignal = ookSignal + ookNoise;

        %demodulation
        ookDemodulated = receivedOOKSignal .* (2 .* carrierSignal);
        ookFiltered = filtfilt(b, a, ookDemodulated);

        %sampling period for demodulation
        samplingPeriod = samplingFreq / dataRate;
        avgPower = amplitude^2/2;
        
        [ookInput, ookOutput] = samplingforthreshold(ookFiltered, samplingPeriod, avgPower, bits);


        ookError = 0;

        for k = 1: bits
            if(ookOutput(k) ~= generatedData(k))
                ookOutput = ookOutput + 1;
            end
        end

        ookError = ookError./bits;
        ookAverageError = ookError + ookAverageError;

    end

    %Plot the 5db SNR signals
    if (snrDB(i) == 5)
    %Plot of original data with respect to time
       figure(2)
       plot(data, 'b');
       title("Original Data")
       xlim([0 1024])

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
    end

end    

% Calculate OOK coherent
e1OOK = (1 / 2) * amplitude^2 / bits;
e0OOk = 0;
ebOOK = (1 / 2) * (e1OOK + e0OOk);
noOOK = ookNoisePower ./ bits ./ 2;
coherentOOK = (1 / 2) .* erfc(sqrt(ebOOK ./ (2 .* noOOK)));
        

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
