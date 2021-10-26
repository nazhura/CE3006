function [input, output] = sampleAndThreshold(filter, period, threshold, bits)
    input = zeros(1, bits);
    output = input;
    avgTime = period / 2;

    for i = 1: bits
        input(i) = filter((2 * (i - 0.5)) * avgTime);
        if(input(i) > threshold)
            output(i) = 1;
        else
            output(i) = 0;
        end
    end
end