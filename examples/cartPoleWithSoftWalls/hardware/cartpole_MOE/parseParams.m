function [layerRanges, layers] = parseParams(layerSize)
    inputNum=5;
    layerNum = length(layerSize);
    layers = cat(2, inputNum, layerSize);
    layerRanges = zeros(layerNum, 4);
    % parse the weights and biases from the vector of parameters
    wbegins = 1;
    for i=2:length(layers)
       wends = wbegins + layers(i-1)*layers(i) - 1;
       bbegins = wends+1;
       bends = wends + layers(i);
       layerRanges(i-1,:) = [wbegins, wends, bbegins, bends];
       wbegins = bends + 1;
    end

end