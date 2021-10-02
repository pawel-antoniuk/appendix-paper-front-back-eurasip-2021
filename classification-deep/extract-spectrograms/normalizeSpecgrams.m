% Paweł Antoniuk 2021
% Bialystok University of Technology

function specgrams = normalizeSpecgrams(specgrams)
    maxVal = max(specgrams, [], [1, 2, 3]);
    minVal = min(specgrams, [], [1, 2, 3]);    
    specgrams = (specgrams - minVal) / (maxVal - minVal);
end