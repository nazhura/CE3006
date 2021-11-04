%% Phase1: Data Generation 0dB to 50dB intervals = 5dB

clear all; close all; clc;

%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 0.2 to 30
SNR_dB = 0:5:50;

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

% Create vector to store all the BER
BER = zeros(1,length(SNR));

% Create vector to store all reconstructed data
reconstructed_data = zeros(length(SNR),bits);

% Set the no of runs of error calculation
noOfRun = 1;
     
% Generate Random Input Binary Data of stated length 1024
% transpose required to convert matrix 1024 x 1 to 1 x 1024
inputBinaryData = transpose(randi([0 1], bits, 1));

% Iterate through all the SNR values
for x = 1:length(SNR)
    sumError = 0;
    for y = 1:noOfRun
        %----Transmission Section----%
        
        % Convert Input Binary Data to +1 and -1
        % If Input Binary Data = 1, 2 * (1 - 0.5) = 1
        % If input Binary Data = 0, 2 * (0 - 0.5) = -1
        inputSignalData = 2 .* (inputBinaryData - 0.5);
        
        % Get Noise Power from SNR
        % SNR = S/N, where assume S = 1
        noisePower = 1 ./ SNR(x);
        
        % Generate noise signal of stated length 1024
        % Noise signal with uniform distribution and generated variance
        generatedNoise = randn(1,bits);
        
        % Noise Signal = sqrt(NoisePower/2) * randn(1,bits)
        noiseSignal = sqrt(noisePower/2) .* generatedNoise;
        
        % Generate Transmission Signal with added Noise
        transmissionSignal = inputSignalData + noiseSignal;
        
        %----Receiver Section----%
        
        % Generate Threshold Logic
        thresholdLevel = 0;
        receivedSignal = zeros(1,bits);
        errorCount = 0;
        
        % Iterate through transmitted signal bits
        for n = 1:bits
            % if transmissionSignal bit > threshold level,
            % append 1 to received signal
            if(transmissionSignal(n) > thresholdLevel)
                receivedSignal(n) = 1;
                % if transmissionSignal bit < threshold level,
                % append 0 to received signal
            elseif(transmissionSignal(n) < thresholdLevel)
                receivedSignal(n) = 0;
            end
            
            % If received bit is not the same as the input bit,
            % increment error count
            if(receivedSignal(n) ~= inputBinaryData(n))
                errorCount = errorCount + 1;
            end
        end
        
        % Compute bit error rate
        errorCount = errorCount ./ bits;
        
        %Accumulate the error of each run
        sumError = errorCount + sumError;
    end
    % Get average error from all the runs and store
    BER(x) = sumError / noOfRun;
    reconstructed_data(x,:) =  receivedSignal;
end

BER_Predicted = 1/2.*erfc(sqrt(SNR));

figure('Name','Measured Data');
semilogy(SNR_dB, BER, 'b-*');
xlabel('SNR_dB');
ylabel('BER');
title('BER against SNR');
hold on
semilogy (SNR_dB,BER_Predicted,'m');
legend('Calculated','Theoretical');
axis([0 50 10^(-10) 1]);
hold off

figure('Name', 'Binary Data Comparison');
subplot(5,1,1);
plot(inputBinaryData);
title("Generated Input Data");
xlim([0 bits])
grid on

subplot(5,1,2);
plot(reconstructed_data(1,:));
title("Reconstructed Data, SNR = 0dB");
xlim([0 bits])
grid on

subplot(5,1,3);
plot(reconstructed_data(2,:));
title("Reconstructed Data, SNR = 5dB");
xlim([0 bits])
grid on

subplot(5,1,4);
plot(reconstructed_data(3,:));
title("Reconstructed Data, SNR = 10dB");
xlim([0 bits])
grid on

subplot(5,1,5);
plot(reconstructed_data(4,:));
title("Reconstructed Data, SNR = 15dB");
xlim([0 bits])
grid on

%% Phase1: Data Generation 0dB to 20dB intervals = 0.2dB

clear all; close all; clc;

%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 0.2 to 30
SNR_dB = 0:0.2:20;

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

% Create vector to store all the BER
BER = zeros(1,length(SNR));

% Create vector to store all reconstructed data
reconstructed_data = zeros(length(SNR),bits);

% Set the no of runs of error calculation
noOfRun = 1;
     
% Generate Random Input Binary Data of stated length 1024
% transpose required to convert matrix 1024 x 1 to 1 x 1024
inputBinaryData = transpose(randi([0 1], bits, 1));

% Iterate through all the SNR values
for x = 1:length(SNR)
    sumError = 0;
    for y = 1:noOfRun
        %----Transmission Section----%
        
        % Convert Input Binary Data to +1 and -1
        % If Input Binary Data = 1, 2 * (1 - 0.5) = 1
        % If input Binary Data = 0, 2 * (0 - 0.5) = -1
        inputSignalData = 2 .* (inputBinaryData - 0.5);
        
        % Get Noise Power from SNR
        % SNR = S/N, where assume S = 1
        noisePower = 1 ./ SNR(x);
        
        % Generate noise signal of stated length 1024
        % Noise signal with uniform distribution and generated variance
        generatedNoise = randn(1,bits);
        
        % Noise Signal = sqrt(NoisePower/2) * randn(1,bits)
        noiseSignal = sqrt(noisePower/2) .* generatedNoise;
        
        % Generate Transmission Signal with added Noise
        transmissionSignal = inputSignalData + noiseSignal;
        
        %----Receiver Section----%
        
        % Generate Threshold Logic
        thresholdLevel = 0;
        receivedSignal = zeros(1,bits);
        errorCount = 0;
        
        % Iterate through transmitted signal bits
        for n = 1:bits
            % if transmissionSignal bit > threshold level,
            % append 1 to received signal
            if(transmissionSignal(n) > thresholdLevel)
                receivedSignal(n) = 1;
                % if transmissionSignal bit < threshold level,
                % append 0 to received signal
            elseif(transmissionSignal(n) < thresholdLevel)
                receivedSignal(n) = 0;
            end
            
            % If received bit is not the same as the input bit,
            % increment error count
            if(receivedSignal(n) ~= inputBinaryData(n))
                errorCount = errorCount + 1;
            end
        end
        
        % Compute bit error rate
        errorCount = errorCount ./ bits;
        
        %Accumulate the error of each run
        sumError = errorCount + sumError;
    end
    % Get average error from all the runs and store
    BER(x) = sumError / noOfRun;
    reconstructed_data(x,:) =  receivedSignal;
end

BER_Predicted = 1/2.*erfc(sqrt(SNR));

figure('Name','Measured Data');
semilogy(SNR_dB, BER, 'b-*');
xlabel('SNR_dB');
ylabel('BER');
title('BER against SNR');
hold on
semilogy (SNR_dB,BER_Predicted,'m');
legend('Calculated','Theoretical');
axis([0 20 10^(-20) 1]);
hold off

figure('Name', 'Binary Data Comparison');
subplot(5,1,1);
plot(inputBinaryData);
title("Generated Input Data");
xlim([0 bits])
grid on

subplot(5,1,2);
plot(reconstructed_data(1,:));
title("Reconstructed Data, SNR = 0dB");
xlim([0 bits])
grid on

subplot(5,1,3);
plot(reconstructed_data(2,:));
title("Reconstructed Data, SNR = 5dB");
xlim([0 bits])
grid on

subplot(5,1,4);
plot(reconstructed_data(3,:));
title("Reconstructed Data, SNR = 10dB");
xlim([0 bits])
grid on

subplot(5,1,5);
plot(reconstructed_data(4,:));
title("Reconstructed Data, SNR = 15dB");
xlim([0 bits])
grid on
