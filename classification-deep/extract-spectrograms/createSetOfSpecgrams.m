% Paweł Antoniuk 2021
% Bialystok University of Technology

function [timeArray, freqArray, specgrams] = ...
        createSetOfSpecgrams(stereoSignal, params) 
	
	stereoSignal = stereoSignal - ones(size(stereoSignal)) * diag(mean(stereoSignal)); % Remove DC offset  
    
    leftChannelSignal = stereoSignal(:, 1);
    rightChannelSignal = stereoSignal(:, 2);
    midSignal = stereoSignal * [1; 1];
    diffSignal = stereoSignal * [1; -1];
    
    midSignal = midSignal - mean(midSignal);
    diffSignal = diffSignal - mean(diffSignal);
    
    specgrams = zeros(params.NTimeFrames, params.NChannels, 4);
    
    [timeArray, freqArray, specgrams(:, :, 1)] = createSpecgram(...
        leftChannelSignal, params);
    [~, ~, specgrams(:, :, 2)] = createSpecgram(rightChannelSignal, ...
        params);
    [~, ~, specgrams(:, :, 3)] = createSpecgram(midSignal, params);
    [~, ~, specgrams(:, :, 4)] = createSpecgram(diffSignal, params);
end
