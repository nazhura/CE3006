function receivedSignal = receivedSignal (signal, signalPower, generatedSignal, SNRVal)
    noisePower = signalPower ./ SNRVal;
    noise = sqrt(noisePower) .* generatedSignal;
    receivedSignal = signal + noise;
end