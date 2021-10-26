function ookbpskFiltered= ookbpskDemo(receivedSignal,carrierSignal, b, a)
    nyquistSampling = 2 .* carrierSignal;
    demodulated = receivedSignal .* nyquistSampling;
    ookbpskFiltered = filtfilt(b, a, demodulated);
end