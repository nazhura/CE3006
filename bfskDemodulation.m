function bfskDiff = bfskDemodulation(receivedBFSKSignal,carrierSignalBFSK1, carrierSignalBFSK2, b, a)
    bfskDemodulated1 = receivedBFSKSignal .* (2 .* carrierSignalBFSK1);
    bfskDemodulated2 = receivedBFSKSignal .* (2 .* carrierSignalBFSK2);
    bfskFiltered1 = filtfilt(b, a, bfskDemodulated1);
    bfskFiltered2 = filtfilt(b, a, bfskDemodulated2);
    bfskDiff = bfskFiltered1 - bfskFiltered2;
end