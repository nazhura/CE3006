%% Phase1: Data Generation

clear all; close all; clc;

%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 0.2 to 30
SNR_dB = 0:0.2:30;

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

% Create vector to store all the BER
BER = [];

% Set the no of runs of error calculation
noOfRun = 30;

% Iterate through all the SNR values
for x = 1:length(SNR)
    sumError = 0;
    for y = 1:noOfRun
        %----Transmission Section----%
        
        % Generate Random Input Binary Data of stated length 1024
        % transpose required to convert matrix 1024 x 1 to 1 x 1024
        inputBinaryData = transpose(randi([0 1], bits, 1));
        
        % Convert Input Binary Data to +1 and -1
        % If Input Binary Data = 1, 2 * 1 - 1 = 1
        % If input Binary Data = 0, 2 * 0 - 1 = -1
        inputSignalData = 2 .* inputBinaryData - 1;
        
        % Get Noise Power from SNR
        % SNR = S/N, where assume S = 1
        noisePower = 1 ./ SNR(x);
        
        % Generate noise signal of stated length 1024
        % Noise signal with uniform distribution and generated variance
        % Noise Signal = sqrt(NoisePower/2) * randn(1,bits)
        noiseSignal = sqrt(noisePower/2) .* randn(1,bits);
        
        % Generate Transmission Signal with added Noise
        transmissionSignal = inputSignalData + noiseSignal;
        
        %----Receiver Section----%
        
        % Generate Threshold Logic
        thresholdLevel = 0;
        receivedSignal = [];
        errorCount = 0;
        
        % Iterate through transmitted signal bits
        for n = 1:bits
            % if transmissionSignal bit > threshold level,
            % append 1 to received signal
            if(transmissionSignal(n) > thresholdLevel)
                receivedSignal = [receivedSignal, 1];
                % if transmissionSignal bit < threshold level,
                % append 0 to received signal
            elseif(transmissionSignal(n) < thresholdLevel)
                receivedSignal = [receivedSignal, 0];
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
end

BER_Predicted = 1/2.*erfc(sqrt(SNR));

figure('Name','Measured Data');
semilogy(SNR_dB, BER, 'b*');
xlabel('SNR_dB');
ylabel('BER');
title('BER against SNR');
hold on
semilogy (SNR_dB,BER_Predicted,'m');
legend('Simulation','Prediction');
axis([0 30 10^(-5) 1]);
hold off

%% Phase1: Data Generation without multiple runs

clear all; close all; clc;

%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 5
SNR_dB = 0:0.01:20;

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

BER = [];

% Iterate through all the SNR values
for x = 1:length(SNR)
    %----Transmission Section----%
    
    % Generate Random Input Binary Data of stated length 1024
    % transpose required to convert matrix 1024 x 1 to 1 x 1024
    inputBinaryData = transpose(randi([0 1], bits, 1));
    
    % Convert Input Binary Data to +1 and -1
    % If Input Binary Data = 1, 2 * 1 - 1 = 1
    % If input Binary Data = 0, 2 * 0 - 1 = -1
    inputSignalData = 2 .* inputBinaryData - 1;
    
    % Get Noise Power from SNR
    % SNR = S/N, where assume S = 1
    noisePower = 1 ./ SNR(x);
    
    % Generate noise signal of stated length 1024
    % Noise signal with uniform distribution and generated variance
    % Noise Signal = sqrt(NoisePower/2) * randn(1,bits)
    noiseSignal = sqrt(noisePower/2) .* randn(1,bits);
    
    % Generate Transmission Signal with added Noise
    transmissionSignal = inputSignalData + noiseSignal;
    
    %----Receiver Section----%
    
    % Generate Threshold Logic
    thresholdLevel = 0;
    receivedSignal = [];
    errorCount = 0;
    
    for n = 1:bits
        % if transmissionSignal bit > threshold level,
        % append 1 to received signal
        if(transmissionSignal(n) > thresholdLevel)
            receivedSignal = [receivedSignal, 1];
            % if transmissionSignal bit < threshold level,
            % append 0 to received signal
        elseif(transmissionSignal(n) < thresholdLevel)
            receivedSignal = [receivedSignal, 0];
        end
        
        % If received bit is not the same as the input bit,
        % increment error count
        if(receivedSignal(n) ~= inputBinaryData(n))
            errorCount = errorCount + 1;
        end
    end
    
    % Compute bit error rate
    BER(x) = errorCount ./ bits;
end

BER_Predicted = 1/2.*erfc(sqrt(SNR));

figure('Name','Measured Data');
semilogy(SNR_dB, BER, 'b*');
xlabel('SNR_dB');
ylabel('BER');
title('BER against SNR');
hold on
semilogy (SNR_dB,BER_Predicted,'m');
legend('Simulation','Prediction');
axis([0 20 10^(-5) 1]);
hold off
%%
clear all; close all;
%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 5
SNR_dB = (0:5:50);

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

%----Transmission Section----%

% Generate Random Input Binary Data of stated length 1024
% transpose required to convert matrix 1024 x 1 to 1 x 1024
inputBinaryData = transpose(randi([0 1], bits, 1));

% Convert Input Binary Data to +1 and -1
% If Input Binary Data = 1, 2 * 1 - 1 = 1
% If input Binary Data = 0, 2 * 0 - 1 = -1
inputSignalData = 2 .* inputBinaryData - 1;

% Get Noise Power from SNR
% SNR = S/N, where assume S = 1
noisePower = 1 ./ SNR(1);

% Generate noise signal of stated length 1024
% Noise signal with uniform distribution and generated variance
% Noise Signal = sqrt(NoisePower/2) * randn(1,bits)
noiseSignal = sqrt(noisePower/2) .* randn(1,bits);

% Generate Transmission Signal with added Noise
transmissionSignal = inputSignalData + noiseSignal;

%%
clear all; close all;
receive = [-1 1 -1 -1 1 1 -1 1 1 -1];
input = [1 0 0 0 0 0 0 0 0 0];
error = 0;
receivedSignal = [];
thresholdLevel = 0;

for n = 1:10
    if(receive(n) > thresholdLevel)
        receivedSignal = [receivedSignal, 1];
    elseif(receive(n) < thresholdLevel)
        receivedSignal = [receivedSignal, 0];        
    end
   if(receivedSignal(n) ~= input(n))
       error = error + 1;
   end
end
BER = error / 10;
%%
close all;clear all;
%----Set Variables----%
% Set number of bits for transmission
bits = 1024;

% Set range of SNR_dB values, increment by 5
SNR_dB = (0:5:50);

% Get SNR Values from SNR_dB
% SNR(dB) = 10log(SNR)
SNR = 10.^(SNR_dB / 10);

% for x = 1:length(SNR)
%     
% end

