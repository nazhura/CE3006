%Setup

clear all;
close all;
clc;

bitCount = 1024;
sigPower = 1;
SNR_dB = 0:1:20;
SNR = power(10, SNR_dB/10);

repeatRuns = 10;
codeword_length = 7;
message_length = 4;

%Run hamming code

for i=1:length(SNR)
    threshold = 0;
    result = zeros(1,repeatRuns);
    received = zeros(1, repeatRuns);
    
    for j=1:repeatRuns
        data = round(rand(1,bitCount));
        signal = times(2, data-1);

        noisePower = rdivide (sigPower, SNR(i));
        noise = sqrt(noisePower/2);
        noise = times(noise, randn(1,bitCount));

        received = signal + noise;

        %Threshold = 0;
        %Error = 0;
    end

    for k = 1: repeatRuns
        result = zeros(1, bitCount);
        for n = 1: bitCount
            if(receive(n)>threshold)
                result(n) = 1;
            else
                result(n) = 0;
            end
        end
    end
    HammingCode = encode(result, codeword_length, message_length, 'hamming/fmt');
end
