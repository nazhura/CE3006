%
% Data Assumptions
% 
% Amplitude of OOK is set to 8
% Hamming Code: Codeword 6 bits, Data 4 bits
% Sampling Frequency = 16 x Carrier Frequency
% 
% Generate data -> Hamming Code Encoding -> OOK/BPSK Modulation -> Noise -> Demodulation -> Comparison
% 

%carrier frequency:
carrierFreq = 10000; %10kHz
carrierFreqFSK1 = 10000;
carrierFreqFSK2 = 40000;

%Self-defined: Codeword length (n) & Message length (k)
codeword_length = 6;
message_length = 4;

%carrier signal 16 times oversampled:
samplingFreq = 16* carrierFreq; %sampling frequency is 16 times the carrier frequency

%baseband data rate:
dataRate = 1000; %1kbps

%number of databits:
bits = 1024;
msgbits = bits*codework_length/message_length;

%sampling rate = sampling frequency / dataRate:
samplingRate = samplingFreq / dataRate;

%amplitude for
amplitude = 8;

%timescale in seconds for ....
time = msgbits/dataRate; %get the time in seconds
period = 1/samplingFreq;
timeScale = 0 : period : time;

%assume a 6th order Butterworth filter with 0.2 normalised cutoff freq
[b, a] = butter(6, 0.2);    %low-pass filter - set in Phase 2

%Carrier Signal Generation
carrierSignal = amp .* cos(2*pi*carrierFreq*timeScale);
carrierSignalFSK1 = amp .* cos(2*pi*carrierFreqFSK1*timeScale);
carrierSignalFSK2 = amp .* cos(2*pi*carrierFreqFSK2*timeScale);

%signal length is for everywhere
signalLength = samplingFreq*msgbits/dataRate + 1;
actualSignalLength = samplingFreq* bits/dataRate + 1;

%SNR
snrDB = 0:0.2:20;
SNR = (10.^(snrDB/10));

%number of test per samples
testSamples = 100;

ookErrorRate = zeros(length(SNR));
ookErrorRate_org = zeros(length(SNR));
bpskErrorRate = zeros(length(SNR));
bfskErrorRate = zeros(length(SNR));
        
%generate data - Hamming
dataHamming = round(rand(1, bits));

%encoding - hamming
encodedHamming = encode(dataHamming, codeword_length, message_length ,'hamming/binary');

%sampling
dataSignal = zeros(1, signalLength);
actualDataSignal = zeros(1, actualSignalLength);
for k = 1: signalLength - 1
    signal(k) = encodedHamming(ceil(k*dataRate/samplingFreq));
end
for k = 1: actualSignalLength - 1
    actualDataSignal(k) = dataSignal(ceil(k*dataRate/samplingFreq));
end
dataSignal(signalLength) = dataSignal(signalLength - 1);
actualDataSignal(actualSignalLength) = actualDataSignal(actualSignalLength - 1);

%==== OOK ====%
%encoded signal
ookSignal = carrierSignal .* dataSignal;
ookSignalPower = (norm(ookSignal)^2)/signalLength;
ookNoisePower = ookSignalPower ./ SNR;

%unencoded signal
orig_ookSignal = carrierSignal(1:actualSignalLength) .* actualDataSignal;
orig_ookSignalPower = norm(orig_ookSignal)^2)/actualSignalLength;
orig_ookNoisePower = orig_ookSignalPower ./ SNR;

%==== BPSK ====%
bpskSourceSignal = dataSignal .* 2 - 1;
bpskSignal = carrierSignal .& bpskSourceSignal;
bpskSignalPower = (norm(bpskSignal)^2)/signalLength;
bpskNoisePower = bpskSignalPower ./ SNR;

%==== BFSK ====%
bfskSourceSignal_high = carrierSignal_FSK1 .* (signal == 1);
bfskSourceSignal_low = carrierSignal_FSK2 .* (signal == 0);
bfskSignal = bfskSourceSignal_low + bfskSourceSignal_high;
bfskSignalPower = (norm(bfskSignal)^2/signalLength;
bfskNoisePower = bfskSignalPower ./ SNR;

% For different SNR values, test over 20 samples
for i = 1 : length(SNR)
    ookAvgError = 0;
    orig_ookAvgError = 0;
    bpskAvgError = 0;
    bfskAvgError = 0;
~ stopped here ~
    
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
