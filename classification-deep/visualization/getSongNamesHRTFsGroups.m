% Paweł Antoniuk 2021
% Bialystok University of Technology

function [songNames, HRTFs, HRTFGroups] = getSongNamesHRTFsGroups(filenames)
    songNames = cell(1, length(filenames));
    HRTFs = cell(1, length(filenames));
    HRTFGroups = cell(1, length(filenames));
    
    for iFilename = 1:length(filenames)
        filename = filenames(iFilename);
        fParts = strsplit(filename, '_');
        songNames{iFilename} = fParts{1};
        HRTFs{iFilename} = fParts{2};
        HRTFGroups{iFilename} = fParts{3};
    end
end