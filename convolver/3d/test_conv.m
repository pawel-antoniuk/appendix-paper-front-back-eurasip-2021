% Paweł Antoniuk 2021
% Bialystok University of Technology

params.TracksDir = 'D:\Recording\RMSNorm 48kHz\Test\AllTheGinIsGone';
params.HRTF = 'HRTF\hutubs\pp1_HRIRs_measured.sofa';
params.outDir = 'test_out';

trackFilenames = dir(fullfile(params.TracksDir, '*.wav'));
tracks = zeros(length(trackFilenames), 504000);
for iTrack = 1:length(trackFilenames)
    trackFilename = trackFilenames(iTrack);
    trackFullFilename = fullfile(trackFilename.folder, trackFilename.name);
    [y,Fs] = audioread(trackFullFilename);
    tracks(iTrack, :) = y;
end

SOFAstart();
SOFA = SOFAload(params.HRTF);

spatTracks = zeros();
for iTrack = 1:size(tracks, 1)
    spatTrack = SOFAspat(tracks(iTrack, :), SOFA, 0, 0);
end
