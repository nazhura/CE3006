function receivedSignal = receivedSignal (signal, signalPower, generatedNoise, SNRVal)
    noisePower = signalPower ./ SNRVal;
    noise = sqrt(noisePower) .* generatedNoise;
    receivedSignal = signal + noise;
end