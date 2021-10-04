% Sławomir Zieliński 2021
% Bialystok University of Technology

function [fileNames, N] = getFileNames(pth)
    folderContent = dir(pth);
    fileSize = cell2mat({folderContent.bytes});
    fileNames = {folderContent(fileSize > 100).name}';
    N = length(fileNames);
end
