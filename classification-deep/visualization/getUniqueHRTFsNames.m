% Paweł Antoniuk 2021
% Bialystok University of Technology

function uniqueHRTFNames = getUniqueHRTFsNames(trainDsFilenames)
    trainDsFilenamesParts = cellfun(@(x) strsplit(x, '_'), trainDsFilenames, 'UniformOutput', false);
    trainDsFilenamesHRTF = cellfun(@(x) x{2}, trainDsFilenamesParts, 'UniformOutput', false);
    uniqueHRTFNames = unique(trainDsFilenamesHRTF);
end