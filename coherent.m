function coherentDetect = coherent(check, noisePower, amplitude, bits)
    e1 = (1 / 2) * amplitude^2 / bits;
    
    % BFSK
    if(check == 1)
        e0 = e1;
        eb = (e0+e1)/2;
        n0 = noisePower ./ bits ./ 2;
        coherentDetect = 0.5 .* erfc(sqrt(eb ./ (2 .* n0)));
    % OOK  
    elseif(check == 0)
        e0 = 0;
        eb = (e0+e1)/2;
        n0 = noisePower ./ bits ./ 2;
        coherentDetect = 0.5 .* erfc(sqrt(eb ./ (2 .* n0)));
    
    % BPSK
    else
        e0 = e1;
        eb = (e0+e1)/2;
        n0 = noisePower ./ bits ./ 2;
        coherentDetect = 0.5 .* erfc(sqrt(eb ./  n0));
    end
   
end

