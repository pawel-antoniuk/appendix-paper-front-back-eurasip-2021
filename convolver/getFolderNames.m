% Sławomir Zieliński 2021
% Bialystok University of Technology

function [folderNames, M] = getFolderNames(pth)
    dirContent = dir(pth);
    isFolder = [dirContent.isdir];
    folderNames = {dirContent(isFolder).name}';
    folderNames(ismember(folderNames,{'.','..','.AppleDouble'})) = [];
    M = length(folderNames);
end
