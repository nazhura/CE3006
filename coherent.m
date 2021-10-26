function coherentDetect = coherent(check, noisePower, amplitude, bits)
    e1 = (1 / 2) * amplitude^2 / bits;
    if(check == 1)
        e0 = e1;
    else
        e0 = 0;
    end
    eb = (e0+e1)/2;

    n0 = noisePower ./ bits ./ 2;

    coherentDetect = 0.5 .* erfc(sqrt(eb ./ (2 .* n0)));
end
