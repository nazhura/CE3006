%% Reference to TheMarvelousWhale - To update


%--Admin stuff--%
clear all; close all; clc;


%{
        This code aims to encode the original data in two different coding
        methods and compare their coding efficiency with the no encode
        version. 
        Coding method chosen: 
            1> Cyclic
            2> Hamming
%}

%define carrier frequency
fc = 10000; %10kHz
%16 times oversampled -> sample freq = 16 fc
fs = 16 * fc;

%define data rate of 1kbps
dataRate = 1000;
%define number of data bits
nBits = 1024;
Enc_nBits = nBits/4*7;   %we doing (7,4) code
%define sampling rate
samplingPeriod = fs / dataRate;

%define Amplitude
Amp = 5;
%define time steps
t = 0: 1/fs : Enc_nBits/dataRate;
t_pure = 0:1/fs : nBits/dataRate;    %pure -- no encode version

%define 6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b_low,a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');


%generate carrier frequency
Carrier = Amp .* cos(2*pi*fc*t);
Carrier_pure = Amp.* cos(2*pi*fc*t_pure);

%calculate signal length
SignalLength = fs*Enc_nBits/dataRate + 1;
SignalLength_pure = fs*nBits/dataRate +1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = -30:1:10;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
%TO PLOT SIGNAL,NOISE and RECEIVE
plot_SNR_dB = 15;



%set run times
Total_Run = 10;

%define placeholder for error calculation
Error_Rate_Hamming = zeros(length(SNR));
Error_Rate_Cyclic = zeros(length(SNR));
Error_Rate_NoEncode = zeros(length(SNR));

%for each SNR value
for i = 1 : length(SNR)
	Avg_Error_Hamming = 0;
    Avg_Error_Cyclic = 0;
    Avg_Error_NoEncode = 0;
    
    %for eacƒh SNR value, average the error over %Total_Run times
	for j = 1 : Total_Run
        
        %-----Data generation-----%
        Data = round(rand(1,nBits));
        EncodeHamming = encode(Data, 7, 4, 'hamming/fmt'); 
        EncodeCyclic = encode(Data,7,4,'cyclic');

        %fill the data stream
        DataStream_Hamming = zeros(1, SignalLength);
        DataStream_Cyclic = zeros(1, SignalLength);

        for k = 1: SignalLength - 1
            DataStream_Hamming(k) = EncodeHamming(ceil(k*dataRate/fs));
            DataStream_Cyclic(k) = EncodeCyclic(ceil(k*dataRate/fs));

        end
        DataStream_Hamming(SignalLength) = DataStream_Hamming(SignalLength - 1);
        DataStream_Cyclic(SignalLength) = DataStream_Cyclic(SignalLength-1);

        DataStream_NoEncode = zeros(1, SignalLength_pure);
        for k = 1:SignalLength_pure -1
            DataStream_NoEncode(k)= Data(ceil(k*dataRate/fs));
        end
        DataStream_NoEncode(SignalLength_pure) = DataStream_NoEncode(SignalLength_pure-1);
        
        
        %----- OOK -----%
        resultOOK_Hamming = OOK_transmission(DataStream_Hamming,SNR(i),Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp);
        resultOOK_Cyclic = OOK_transmission(DataStream_Cyclic,SNR(i),Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp);
        resultOOK_NoEncode = OOK_transmission(DataStream_NoEncode,SNR(i),Carrier_pure,SignalLength_pure,samplingPeriod,nBits,Amp);
        
        decodedHamming = decode(resultOOK_Hamming,7,4,'hamming/fmt');
        decodedCyclic = decode(resultOOK_Cyclic,7,4,'cyclic');
        
        %--Calculate Error--%
        ErrorHamming = 0;
        ErrorCyclic = 0;
        ErrorNoEncode = 0;
        for k = 1: nBits
            if(decodedHamming(k) ~= Data(k))
                ErrorHamming = ErrorHamming + 1;
            end
            if(decodedCyclic(k) ~= Data(k))
                ErrorCyclic = ErrorCyclic + 1;
            end
            if (resultOOK_NoEncode(k) ~= Data(k));
                ErrorNoEncode = ErrorNoEncode + 1;
            end
        end
        Avg_Error_Hamming = ErrorHamming + Avg_Error_Hamming;
        Avg_Error_Cyclic = ErrorCyclic + Avg_Error_Cyclic;
        Avg_Error_NoEncode = ErrorNoEncode + Avg_Error_NoEncode;
    end
    Error_Rate_Hamming(i) = Avg_Error_Hamming/Total_Run/nBits;
    Error_Rate_Cyclic(i) = Avg_Error_Cyclic/Total_Run/nBits;
    Error_Rate_NoEncode(i) = Avg_Error_NoEncode/Total_Run/nBits;
end


%Error plot
figure(1);
semilogy (SNR_dB, Error_Rate_Hamming,'r-*');
hold on
semilogy(SNR_dB, Error_Rate_Cyclic, 'b-*');
hold off
hold on
semilogy(SNR_dB, Error_Rate_NoEncode, 'k-*');
hold off
title('Error rate of cyclic and hamming and No Encoding for different SNR');
legend('hamming','cyclic','None');
ylabel('Pe');
xlabel('Eb/No')




%%--HELPER FUNCTION--%%
function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end


%This function simulates the decision device
function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end

%This function is a wrapper for Phase 2 OOK 
function result_OOK = OOK_transmission(DataStream,SNR_val,Carrier,SignalLength,samplingPeriod,Enc_nBits,Amp)
        Signal_OOK = Carrier .* DataStream;
        [b_low,a_low] = butter(6, 0.2);
        %generate noise 
        Signal_Power_OOK = (norm(Signal_OOK)^2)/SignalLength;  %Sum of squared signal amp over signal length
		Noise_Power_OOK = Signal_Power_OOK ./SNR_val;
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		
        %transmission
		ReceiveOOK = Signal_OOK+NoiseOOK;
        %detection -- square law device
        SquaredOOK = ReceiveOOK .* ReceiveOOK;
        %low pass filter
        FilteredOOK = filtfilt(b_low, a_low, SquaredOOK);
         
        %sample and decision device
        sampledOOK = sample(FilteredOOK, samplingPeriod, Enc_nBits);
        result_OOK = decision_device(sampledOOK,Enc_nBits, Amp/2);  %--OOK threshold is 0.5*(A+0)

end

%% Reference to tanhauhau - to update

clear all; close all; clc;
%Define Signal length
Num_Bit = 1024;
%Define signal power                   
Signal_Power = 1; 
%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

% Set run times
Total_Run = 20;
%Different SNR value
for i = 1 : length(SNR)
	%Avg_Error = 0;
    %Error_Rate=zeros(1,Total_Run);
    result=zeros(1,Total_Run);
    Receive=zeros(1,Total_Run);
    threshold=0;
    
	for j = 1 : Total_Run
		%Input singal
		%Generate random binary digit(0 or 1)
		Data = round(rand(1,Num_Bit));
		%Convert binary digit to (-1 or +1)
		Signal = 2 .* Data - 1;
		%Generate Noise  		
		Noise_Power = Signal_Power ./SNR(i);
		Noise = sqrt(Noise_Power/2) .*randn(1,Num_Bit);   
		%Received Signal 
		Receive = Signal+Noise; 

		Threshold = 0;
		Error = 0; 
		%for k= 1 : Num_Bit
            
		%	if (Receive(k)>= Threshold) && Data(k)==0||(Receive(k)<Threshold && Data(k)==1)
		%		Error = Error+1;
		%	end
		%end
		%Calculate bit error rate
		%Error = Error ./Num_Bit;  
		%Calculate the average error for every runtime		
		%Avg_Error = Error + Avg_Error;   
        
	end
	%Error_Rate(i) = Avg_Error / Total_Run;
    for k=1: Total_Run
        result=zeros(1,Num_Bit);
        for n= 1: Num_Bit
            if(Receive(n)>threshold)
                result(n)=1;
            else
                result(n)=0;
            end
        end
    end
        %result(k) = HammingEncodingThreshold(Num_Bit, threshold, Receive);
    
    EncodeHamming= encode(result,7,4,'hamming/fmt');
    
end

%% Reference to lolfuljames - cyclic
clear all; close all; clc;

% Codeword length (n) & Message length (k)
codeword_length = 7;
message_length = 4;

% Data length (1024) & Encoded data length (1024/4 * 7 = 1792)
signal_length = 1024;
encoded_signal_length = signal_length/message_length * codeword_length;

% Carrier parameters and signal (amplitude = 5) as required
carrier_freq = 10000;
sampling_freq = 16 * carrier_freq;
data_rate = 1000;
t = 0: 1/sampling_freq : encoded_signal_length/data_rate;
amp = 5;

%For FSK modulation
fsk_freq_1 = 30000;
fsk_freq_2 = 10000;

%Carrier signal generation
carrier_signal = amp .* cos(2*pi*carrier_freq*t);
fsk_carrier_signal_1 = amp .* cos(2*pi*fsk_freq_1*t);
fsk_carrier_signal_2 = amp .* cos(2*pi*fsk_freq_2*t);

% Sampled signal length
sampled_signal_length = sampling_freq*encoded_signal_length/data_rate + 1;
sampled_unencoded_signal_length = sampling_freq*signal_length/data_rate + 1;

% Number of samples
nb_samples = 200;

% Low-pass filter - 6th order, 0.2 cutoff frequency
[b, a] = butter(6, 0.2);

% maximum SNR
MAX_dB = 20;

% Generator and Syndrome Table for Cyclic Code (Channel Encoding-Decoding)
% Generator polynomial & Parity-check matrix for cyclic encoding
genpoly = cyclpoly(codeword_length, message_length);
parmat = cyclgen(codeword_length, genpoly);
% Syndrome table for cyclic decoding
trt = syndtable(parmat);


% SNR values to test
SNR_dB = 0:1:10;
SNR = (10.^(SNR_dB/10));

OOK_error_rate = zeros([length(SNR) 1]);
unencoded_OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);
BFSK_error_rate = zeros([length(SNR) 1]);

% Original signal 
original_signal = round(rand(1,signal_length));
% Cyclic encoding
encoded_signal = encode(original_signal, codeword_length, message_length, 'cyclic/binary', genpoly);

% Sampling
sampled_signal = zeros(1, sampled_signal_length);
sampled_unencoded_signal = zeros(1, sampled_unencoded_signal_length);
for k = 1: sampled_signal_length - 1
    sampled_signal(k) = encoded_signal(ceil(k*data_rate/sampling_freq));
end
for k = 1: sampled_unencoded_signal_length - 1
    sampled_unencoded_signal(k) = original_signal(ceil(k*data_rate/sampling_freq));
end
sampled_signal(sampled_signal_length) = sampled_signal(sampled_signal_length - 1);
sampled_unencoded_signal(sampled_unencoded_signal_length) = sampled_unencoded_signal(sampled_unencoded_signal_length - 1);

% Modulation: On-Off Keying
OOK_signal = carrier_signal .* sampled_signal;
unencoded_OOK_signal = carrier_signal(1:sampled_unencoded_signal_length) .* sampled_unencoded_signal;

% Modulation: Binary Phase Shift Keying 
BPSK_source_signal = sampled_signal .* 2 - 1; % put to -1 +1
BPSK_signal = carrier_signal .* BPSK_source_signal;

%Modulation: Binary FSK
BFSK_source_signal_1 = fsk_carrier_signal_1 .* (sampled_signal == 1);
BFSK_source_signal_0 = fsk_carrier_signal_2 .* (sampled_signal == 0);
BFSK_signal = BFSK_source_signal_1 + BFSK_source_signal_0;

OOK_signal_power = (norm(OOK_signal)^2)/sampled_signal_length;
unencoded_OOK_signal_power = (norm(unencoded_OOK_signal)^2)/sampled_unencoded_signal_length;
BPSK_signal_power = (norm(BPSK_signal)^2)/sampled_signal_length;
BFSK_signal_power = (norm(BFSK_signal)^2)/sampled_signal_length;

% For each value of SNR, test of 20 samples.
for i = 1 : length(SNR) 
    
	OOK_average_error = 0;
    unencoded_OOK_average_error = 0;
	BPSK_average_error = 0;
    BFSK_average_error = 0;
    result = zeros(1, nb_samples);
    
    for j = 1 : nb_samples
        
        % Generate White Gaussian Channel Noise for OOK and BPSK
        noise_power_OOK = OOK_signal_power ./ SNR(i);
        noise_OOK = sqrt(noise_power_OOK) .*randn(1,sampled_signal_length);
        
        noise_power_unencoded_OOK = unencoded_OOK_signal_power ./ SNR(i);
        noise_unencoded_OOK =  sqrt(noise_power_unencoded_OOK) .*randn(1,sampled_unencoded_signal_length);

        noise_power_BPSK = BPSK_signal_power ./ SNR(i);
        noise_BPSK = sqrt(noise_power_BPSK) .*randn(1,sampled_signal_length);
        
        noise_power_BFSK = BFSK_signal_power ./ SNR(i);
        noise_BFSK = sqrt(noise_power_BFSK) .*randn(1,sampled_signal_length);
        
        % Received Signal with added Channel Noise for OOK and BPSK
        OOK_received = OOK_signal + noise_OOK;
        unencoded_OOK_received = unencoded_OOK_signal + noise_unencoded_OOK;
        BPSK_received = BPSK_signal + noise_BPSK;
        BFSK_received = BFSK_signal + noise_BFSK;
        
        % Non-Coherent Detection: OOK (Lecture notes 04 - pg 18)
        % Squared
        OOK_squared = OOK_received .* 2 .* carrier_signal;
        unencoded_OOK_squared = unencoded_OOK_received .* 2 .* carrier_signal(1:sampled_unencoded_signal_length);

        % Low-Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
        unencoded_OOK_filtered = filtfilt(b, a, unencoded_OOK_squared);
  
        % Non-Coherent Detection: BPSK (Lecture notes 04 - pg 50)
        % Squared
        BPSK_squared = BPSK_received .* 2 .* carrier_signal;
        % Low-Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_squared);
        
        %BFSK detection
        BFSK_carrier_1_corr = BFSK_received .* (2 .* fsk_carrier_signal_1);
        BFSK_branch_1_filtered = filtfilt(b, a, BFSK_carrier_1_corr);
        BFSK_carrier_2_corr = BFSK_received .* (2 .* fsk_carrier_signal_2);
        BFSK_branch_2_filtered = filtfilt(b, a, BFSK_carrier_2_corr);
        BFSK_differenced = BFSK_branch_1_filtered - BFSK_branch_2_filtered;
        
        % Demodulation: Sampling and Threshold
        sampling_period = sampling_freq/data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sampling_period, 25/2, encoded_signal_length);
        [unencoded_OOK_sample, unencoded_OOK_result] = sample_and_threshold(unencoded_OOK_filtered, sampling_period, 25/2, signal_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sampling_period, 0, encoded_signal_length);
        [BFSK_sample, BFSK_result] = sample_and_threshold(BFSK_differenced, sampling_period, 0, encoded_signal_length);

        % Cyclic Decoding
        decoded_signal_OOK = decode(OOK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BPSK = decode(BPSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        decoded_signal_BFSK = decode(BFSK_result, codeword_length, message_length, 'cyclic/binary', genpoly, trt);
        
        % Get Bit Errors
        OOK_error = biterr(decoded_signal_OOK, original_signal) ./ signal_length;
        unencoded_OOK_error = biterr(unencoded_OOK_result, original_signal) ./ signal_length;
        BPSK_error = biterr(decoded_signal_BPSK, original_signal) ./ signal_length;
        BFSK_error = biterr(decoded_signal_BFSK, original_signal) ./signal_length;
                
        OOK_average_error = OOK_error + OOK_average_error;
        unencoded_OOK_average_error = unencoded_OOK_error + unencoded_OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
        BFSK_average_error = BFSK_error + BFSK_average_error;
    end
         
    OOK_error_rate(i) = OOK_average_error / nb_samples;
    BPSK_error_rate(i) = BPSK_average_error / nb_samples;
    BFSK_error_rate(i) = BFSK_average_error / nb_samples;
    unencoded_OOK_error_rate(i) = unencoded_OOK_average_error / nb_samples; 
    % Plots for SNR @ 5dB
    if (i == 5)
        figure(2)
        subplot(2, 1, 1);
        plot(original_signal, 'b');
        title("Original Signal")
        xlim([0 2000])

        subplot(2, 1, 2);
        plot(encoded_signal, 'b');
        title("Cyclic Encoded Data")
        xlim([0 2000])

        figure(3)
        subplot(4, 1, 1);
        spectrogram(OOK_signal,'yaxis')
        title("Transmitted OOK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(OOK_received,'yaxis')
        title("Received OOK Modulated Signal")

        subplot(4, 1, 3);
        plot(OOK_sample)
        title("OOK Demodulated Signal")

        subplot(4, 1, 4);
        plot(decoded_signal_OOK)
        title("OOK Decoded Signal");

        figure(3)
        subplot(4, 1, 1);
        spectrogram(BPSK_signal,'yaxis')
        title("Transmitted BPSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BPSK_received,'yaxis')
        title("Received BPSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BPSK_sample)
        title("BPSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(decoded_signal_BPSK)
        title("BPSK Decoded Signal");
        
        figure(5)
        subplot(4, 1, 1);
        spectrogram(BFSK_signal,'yaxis')
        title("Transmitted BFSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BFSK_received,'yaxis')
        title("Received BFSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BFSK_sample)
        title("BFSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(decoded_signal_BFSK)
        title("BFSK Decoded Signal");

        figure(6);
        subplot(4, 1, 1);
        plot(original_signal);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(4, 1, 2);
        plot(decoded_signal_OOK);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(4, 1, 3);
        plot(decoded_signal_BPSK);
        title("BPSK Decoded Data");
        xlim([0 1024]);
        
        subplot(4, 1, 4);
        plot(decoded_signal_BFSK);
        title("BFSK Decoded Data");
        xlim([0 1024]);
    end

end

% Plot OOK vs BPSK bit error rate
figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
p3 = semilogy(SNR_dB, BFSK_error_rate, 'g-*');
p4 = semilogy(SNR_dB, unencoded_OOK_error_rate, 'k-*');
hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1) p3(1) p4(1)],{'Cyclic/OOK','Cyclic/BPSK',' Cyclic/BFSK','Unencoded/OOK'})
xlim([0 50]);


%% Reference to lolfuljames - hamming.m

%
% Data Assumptions
% 
% Amplitude of OOK is set to 5
% Hamming Code: Codewode 7 bits, Data 4 bits
% Sampling Frequency = 16 x Carrier Frequency
% 
% Generate data -> Hamming Code Encoding -> OOK/BPSK Modulation -> Noise -> Demodulation -> Comparison
% 

carrier_freq = 10000; %10kHz
sample_freq = 16 * carrier_freq;
data_rate = 1000; %1kbps
data_length = 1024;
encoded_signal_length = 1792;
amp = 5;

%For FSK modulation
fsk_freq_1 = 30000;
fsk_freq_2 = 10000;

% Low Pass 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);

% Time simulation
t = 0: 1/sample_freq : encoded_signal_length/data_rate;

% Carrier Signal Generation
carrier_signal = amp .* cos(2*pi*carrier_freq*t);
fsk_carrier_signal_1 = amp .* cos(2*pi*fsk_freq_1*t);
fsk_carrier_signal_2 = amp .* cos(2*pi*fsk_freq_2*t);

% Length of transmitted signal
signal_length = sample_freq*encoded_signal_length/data_rate + 1;
signal_length_unencoded = sample_freq*data_length/data_rate + 1;

% SNR values to test
SNR_dB = 0:1:20;
SNR = (10.^(SNR_dB/10));

% Number of tests per SNR
test_samples = 200;

OOK_error_rate = zeros([length(SNR) 1]);
unencoded_OOK_error_rate = zeros([length(SNR) 1]);
BPSK_error_rate = zeros([length(SNR) 1]);
BFSK_error_rate = zeros([length(SNR) 1]);
        
% Generate Hamming encoded signals
data = round(rand(1,data_length));
hamming_signal= encode(data,7,4,'hamming/binary');

signal = zeros(1, signal_length);
signal_unencoded = zeros(1, signal_length_unencoded);
for k = 1: signal_length - 1
    signal(k) = hamming_signal(ceil(k*data_rate/sample_freq));
end
for k = 1: signal_length_unencoded - 1
    signal_unencoded(k) = data(ceil(k*data_rate/sample_freq));
end
signal(signal_length) = signal(signal_length - 1);
signal_unencoded(signal_length_unencoded) = signal_unencoded(signal_length_unencoded - 1);

% OOK Modulation
OOK_signal = carrier_signal .* signal;
unencoded_OOK_signal = carrier_signal(1:length(signal_unencoded)) .* signal_unencoded;

% BPSK Modulation
BPSK_source_signal = signal .* 2 - 1;
BPSK_signal = carrier_signal .* BPSK_source_signal;

% BFSK Modulation
BFSK_source_signal_1 = fsk_carrier_signal_1 .* (signal == 1);
BFSK_source_signal_0 = fsk_carrier_signal_2 .* (signal == 0);
BFSK_signal = BFSK_source_signal_1 + BFSK_source_signal_0;

OOK_signal_power = (norm(OOK_signal)^2)/signal_length;
unencoded_OOK_signal_power = (norm(unencoded_OOK_signal)^2)/signal_length_unencoded;
BPSK_signal_power = (norm(BPSK_signal)^2)/signal_length;
BFSK_signal_power = (norm(BFSK_signal)^2)/signal_length;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
	OOK_average_error = 0;
    unencoded_OOK_average_error = 0;
    BPSK_average_error = 0;
    BFSK_average_error = 0;
    
	for j = 1 : test_samples
        
        % Generate noise
		noise_power_OOK = OOK_signal_power ./SNR(i);
		noise_OOK = sqrt(noise_power_OOK) .*randn(1,signal_length);
        
        noise_power_unencoded_OOK = unencoded_OOK_signal_power ./ SNR(i);
        noise_unencoded_OOK =  sqrt(noise_power_unencoded_OOK) .*randn(1,signal_length_unencoded);
        
        noise_power_BPSK = BPSK_signal_power ./SNR(i);
        noise_BPSK = sqrt(noise_power_OOK) .*randn(1,signal_length);
        
        noise_power_BFSK = BFSK_signal_power ./SNR(i);
        noise_BFSK = sqrt(noise_power_BFSK) .*randn(1,signal_length);
        
		% OOK Signal on Receiver's end
		OOK_received = OOK_signal+noise_OOK;
        unencoded_OOK_received = unencoded_OOK_signal + noise_unencoded_OOK;
        
        % Start of OOK Detection
        OOK_squared = OOK_received .* 2.* carrier_signal;
        unencoded_OOK_squared = unencoded_OOK_received .* 2 .* carrier_signal(1:length(signal_unencoded));
        
        % Low Pass Filter
        OOK_filtered = filtfilt(b, a, OOK_squared);
        unencoded_OOK_filtered = filtfilt(b, a, unencoded_OOK_squared);
        
		% BPSK Signal on Receiver's end
		BPSK_received = BPSK_signal+noise_BPSK;
        
        % Non-coherent detection
        BPSK_squared = BPSK_received .* (2.* carrier_signal);
        % Low Pass Filter
        BPSK_output = filtfilt(b, a, BPSK_squared);
        
        %Received Signal BFSK
        BFSK_received = BFSK_signal+noise_BFSK;
        
        %BFSK detection
        BFSK_carrier_1_corr = BFSK_received .* (2 .* fsk_carrier_signal_1);
        BFSK_branch_1_filtered = filtfilt(b, a, BFSK_carrier_1_corr);
        BFSK_carrier_2_corr = BFSK_received .* (2 .* fsk_carrier_signal_2);
        BFSK_branch_2_filtered = filtfilt(b, a, BFSK_carrier_2_corr);
        BFSK_differenced = BFSK_branch_1_filtered - BFSK_branch_2_filtered;
        
        % Demodulation by sample & threshold
        sample_period = sample_freq / data_rate;
        [OOK_sample, OOK_result] = sample_and_threshold(OOK_filtered, sample_period, (amp^2)/2, encoded_signal_length);
        [unencoded_OOK_sample, unencoded_OOK_result] = sample_and_threshold(unencoded_OOK_filtered, sample_period, (amp^2)/2, data_length);
        [BPSK_sample, BPSK_result] = sample_and_threshold(BPSK_output, sample_period, 0, encoded_signal_length);
        [BFSK_sample, BFSK_result] = sample_and_threshold(BFSK_differenced, sample_period, 0, encoded_signal_length);

        OOK_decoded = decode(OOK_result,7,4,'hamming/binary');
        BPSK_decoded = decode(BPSK_result,7,4,'hamming/binary');
        BFSK_decoded = decode(BFSK_result,7,4,'hamming/binary');
        
        OOK_error =  biterr(OOK_decoded, data) ./encoded_signal_length;
        unencoded_OOK_error = biterr(unencoded_OOK_result, data) ./data_length;
        BPSK_error = biterr(BPSK_decoded, data) ./encoded_signal_length;
        BFSK_error = biterr(BFSK_decoded, data) ./encoded_signal_length;
        OOK_average_error = OOK_error + OOK_average_error;
        unencoded_OOK_average_error = unencoded_OOK_error + unencoded_OOK_average_error;
        BPSK_average_error = BPSK_error + BPSK_average_error;
        BFSK_average_error = BFSK_error + BFSK_average_error;
    end
	OOK_error_rate(i) = OOK_average_error / test_samples;
    unencoded_OOK_error_rate(i) = unencoded_OOK_average_error / test_samples;
    BPSK_error_rate(i) = BPSK_average_error / test_samples;
    BFSK_error_rate(i) = BFSK_average_error / test_samples;
    
%     Plot the 5db SNR signals
    if (SNR_dB(i) == 5)
        figure(2)
        subplot(2, 1, 1);
        plot(data, 'b');
        title("Original Data")
        xlim([0 2000])

        subplot(2, 1, 2);
        plot(hamming_signal, 'b');
        title("Hamming Encoded Data")
        xlim([0 2000])

        figure(3)
        subplot(4, 1, 1);
        spectrogram(OOK_signal,'yaxis')
        title("Transmitted OOK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(OOK_received,'yaxis')
        title("Received OOK Modulated Signal")

        subplot(4, 1, 3);
        plot(OOK_sample)
        title("OOK Demodulated Signal")

        subplot(4, 1, 4);
        plot(OOK_decoded)
        title("OOK Decoded Signal");

        figure(4)
        subplot(4, 1, 1);
        spectrogram(BPSK_signal,'yaxis')
        title("Transmitted BPSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BPSK_received,'yaxis')
        title("Received BPSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BPSK_sample)
        title("BPSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(BPSK_decoded)
        title("BPSK Decoded Signal");
        
        figure(5)
        subplot(4, 1, 1);
        spectrogram(BFSK_signal,'yaxis')
        title("Transmitted BFSK Modulated Signal")

        subplot(4, 1, 2);
        spectrogram(BFSK_received,'yaxis')
        title("Received BFSK Modulated Signal")

        subplot(4, 1, 3);
        plot(BFSK_sample)
        title("BFSK Demodulated Signal")

        subplot(4, 1, 4);
        plot(BFSK_decoded)
        title("BFSK Decoded Signal");

        figure(6);
        subplot(4, 1, 1);
        plot(data);
        title("Original Data");
        xlim([0 1024]);
        ylim([0 1]);

        subplot(4, 1, 2);
        plot(OOK_decoded);
        title("OOK Decoded Data");
        xlim([0 1024]);

        subplot(4, 1, 3);
        plot(BPSK_decoded);
        title("BPSK Decoded Data");
        xlim([0 1024]);
        
        subplot(4, 1, 4);
        plot(BFSK_decoded);
        title("BFSK Decoded Data");
        xlim([0 1024]);
    end
end

% Plot OOK vs DBSK bit error rate
figure(1)
p1 = semilogy(SNR_dB, OOK_error_rate,'r-*');
hold on
p2 = semilogy(SNR_dB, BPSK_error_rate, 'b-*');
p3 = semilogy(SNR_dB, BFSK_error_rate, 'g-*');
p4 = semilogy(SNR_dB, unencoded_OOK_error_rate, 'k-*');

hold off
ylabel('Bit Error Rate (BER)');
xlabel('SNR (dB)');
legend([p1(1) p2(1) p3(1) p4(1)],{'Hamming/OOK','Hamming/BPSK',' Hamming/BFSK','Unencoded/OOK'})
xlim([0 50]);

