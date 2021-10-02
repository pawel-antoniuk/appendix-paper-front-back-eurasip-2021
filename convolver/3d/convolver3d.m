% Paweł Antoniuk 2021
% Bialystok University of Technology

%% Initialize
clearvars; close all; clc;
SOFAstart;
SOFAgetVersion();

%% Params
params.HRTFBaseDir = 'HRTF';
params.RecordingsBaseDir = 'recordings\';
params.RecordingsScenarios = ["Train", "Test"];
params.SpatOutputDir = 'binaural\spat';
params.MetaOutputDir = 'binaural\meta';
params.FinalResultsOutputDir = 'binaural';
params.RecordingsExpectedFs = 48000;
params.RecordingLoadRange = [2.5 inf];
params.RecordingSpatRange = [0.5 7];
params.RecordingFadeTime = [0.01 0.01];
params.RecordingLevelScale = 0.9;
params.NChannels = 2;
params.EnsambleAngleWidthEps = 1e-4;
params.EnsambleAngleWidths = [15 30 45 60 90];
params.EnsambleDirections = [
    0 0
    -180 0
    0 90
    0 -90];
params.EnsambleSceneNames = [
    "front"
    "back"
    "up"
    "down"];

%% Load HRTFs and get audio filenames
HRTFs = loadHRTFs(params);
audioFilenames = getAudioFilenames(params);

%% Spatialize songs
tic
nAudioFilenames = length(audioFilenames);
allSpatMetaresults = cell(nAudioFilenames, 1);
parfor iAudioFilename = 1:nAudioFilenames
    %% Load audio tracks and spatialize them
    audioFilename = audioFilenames(iAudioFilename);
    tracks = loadAudioTracks(audioFilename, params);
    [spatResults,spatMetaresults] = spatializeSong(HRTFs, tracks, params);
    allSpatMetaresults{iAudioFilename} = spatMetaresults;
    
    %% Postprocess results
    spatResults = posprocessSpatResults(spatResults, params);
    
    %% Save results
    saveSpatResults(spatResults, spatMetaresults, HRTFs, audioFilename, params);
    
    fprintf("Progress  [audio: %d/%d] (%s)\n", ...
        iAudioFilename, nAudioFilenames, audioFilename.name);
end
toc

%% Plot scenes summary
plotAudioScene(HRTFs, ...
    cell2mat(reshape(allSpatMetaresults, 1, 1, 1, [])), params);

%% Save workspace
save(fullfile(params.FinalResultsOutputDir, 'workspace'), '-v7.3');

%% --- ROUTINES ---
% Load HRTFs routine
function HRTFs = loadHRTFs(params)

HRTFFilenames = dir(fullfile(params.HRTFBaseDir, '*', '*.sofa'));
[~, HRTFIdx] = sort({HRTFFilenames.name});
HRTFFilenames = HRTFFilenames(HRTFIdx);

% HRTF struct definition
HRTFs = struct('Id', [], ...
    'Name', [], ...
    'Folder', [], ...
    'HRTFGroup', [], ...
    'SOFA', [], ...
    'Position', [], ...
    'Distance', []);
HRTFGroupData = containers.Map;

for iHRTF = 1:length(HRTFFilenames)
    filename = HRTFFilenames(iHRTF);
    fullFilename = fullfile(filename.folder, filename.name);
    
    HRTFs(iHRTF) = loadHRTF(iHRTF, fullFilename);   
    
    if HRTFs(iHRTF).SOFA.Data.SamplingRate ~= params.RecordingsExpectedFs
        [loadStatus,HRTFs(iHRTF)] = tryLoadResampledHRTF(iHRTF, ...
            HRTFs(iHRTF), params);
        if ~loadStatus
            resampleAndSave(HRTFs(iHRTF), params);
            [loadStatus,HRTFs(iHRTF)] = tryLoadResampledHRTF(iHRTF, ...
                HRTFs(iHRTF), params);
            
            if ~loadStatus
                error('Cannot find previously resampled HRTF');
            end
        end
    end
    
    if ~isKey(HRTFGroupData, HRTFs(iHRTF).HRTFGroup)
        HRTFGroupData(HRTFs(iHRTF).HRTFGroup) = [];
    end
    
    for jHRTF = HRTFGroupData(HRTFs(iHRTF).HRTFGroup)
        if length(HRTFs(iHRTF).Position) ~= length(HRTFs(jHRTF).Position) ...
                || ~all(HRTFs(iHRTF).Position == HRTFs(jHRTF).Position, 'all')
            warning('[%s][%s] Inconsistent source positions with %s', ...
                HRTFs(iHRTF).HRTFGroup, ...
                HRTFs(iHRTF).Name, ...
                HRTFs(jHRTF).Name);
        end
    end

    HRTFGroupData(HRTFs(iHRTF).HRTFGroup) = [...
        HRTFGroupData(HRTFs(iHRTF).HRTFGroup) iHRTF];   
    
    fprintf('[%s][%s] azimuth: [%d, %d]; elevation: [%d, %d]; distance: %d\n', ...
        HRTFs(iHRTF).HRTFGroup, ...
        HRTFs(iHRTF).Name, ...
        min(HRTFs(iHRTF).Position(:, 1)), ...
        max(HRTFs(iHRTF).Position(:, 1)), ...
        min(HRTFs(iHRTF).Position(:, 2)), ...
        max(HRTFs(iHRTF).Position(:, 2)), ...
        HRTFs(iHRTF).Distance);
    
    if length(HRTFs(iHRTF).Distance) > 1
        error('Multiple distances in single HRTF are not supported');
    end
    
    if HRTFs(iHRTF).SOFA.Data.SamplingRate ~= params.RecordingsExpectedFs
        error('[%s][%s] Resampling from %d Hz to %d Hz', ...
            HRTF.HRTFGroup, HRTF.Name, ...
            HRTF.SOFA.Data.SamplingRate, ...
            params.RecordingsExpectedFs);
    end
end

end


% Try load resampled HRTF routine
function [loadStatus,HRTF] = tryLoadResampledHRTF(id, HRTF, params)

resampledSOFAdir = fullfile(params.HRTFBaseDir, ...
    ['_resampled_' num2str(params.RecordingsExpectedFs)], ...
    HRTF.HRTFGroup);
resampledSOFAfilename = ['_resampled_' ...
    num2str(params.RecordingsExpectedFs) '_' HRTF.Name];
fullSOFAfilename = fullfile(resampledSOFAdir, resampledSOFAfilename);

if ~exist(fullSOFAfilename, 'file')
    loadStatus = false;
else
    loadStatus = true;
    HRTF = loadHRTF(id, fullSOFAfilename);
end

end


% Load HRTF routine
function HRTF = loadHRTF(id, filename)

listing = dir(filename);
fullFilename = fullfile(listing.folder, listing.name);
filenameParts = split(listing.folder, filesep);
SOFA = SOFAload(fullFilename);
APV = SOFAcalculateAPV(SOFA);

HRTF.Id = id;
HRTF.Name = listing.name;
HRTF.Folder = listing.folder;
HRTF.HRTFGroup = filenameParts{end};
HRTF.SOFA = SOFA;
% HRTF.Position = HRTF.SOFA.SourcePosition(:, 1:2);
HRTF.Position = APV(:, 1:2);
HRTF.Distance = unique(HRTF.SOFA.SourcePosition(:, 3));

% If the number of samples is odd, remove the last sample
% It fixes the problem with SOFA-based convolution
% (S. Zieliński)
if mod(HRTF.SOFA.API.N, 2) ~= 0
    tmpIR = HRTF.SOFA.Data.IR(:, :, 1:end-1); % Remove last sample
    HRTF.SOFA.Data.IR = tmpIR;
    HRTF.SOFA.API.N = size(tmpIR, 3);
end        

end


% Resample and save routine
function HRTF = resampleAndSave(HRTF, params)

fprintf('[%s][%s] Resampling from %d Hz to %d Hz\n', ...
    HRTF.HRTFGroup, HRTF.Name, ...
    HRTF.SOFA.Data.SamplingRate, ...
    params.RecordingsExpectedFs);

HRTF.SOFA = SOFAresample(HRTF.SOFA, params.RecordingsExpectedFs);

resampledSOFAdir = fullfile(params.HRTFBaseDir, ...
    ['_resampled_' num2str(params.RecordingsExpectedFs)], ...
    HRTF.HRTFGroup);
resampledSOFAfilename = ['_resampled_' ...
    num2str(params.RecordingsExpectedFs) '_' HRTF.Name];

if ~exist(resampledSOFAdir, 'dir')
    mkdir(resampledSOFAdir);
end

fullSOFAfilename = fullfile(resampledSOFAdir, resampledSOFAfilename);
HRTF.SOFA = SOFAsave(fullSOFAfilename, HRTF.SOFA, 0);

end


% Resample SOFA routine
function Obj = SOFAresample(Obj, targetFs)

currentFs = Obj.Data.SamplingRate;

if currentFs == targetFs
    return
end

% Based on HRTFsamplingRateConverter10.m (S. Zieliński)
M = size(Obj.Data.IR,1); % Number of measurements
N = size(Obj.Data.IR,3); % Length of measurements
IR = Obj.Data.IR;
IR2 = zeros(M, 2, round(targetFs / currentFs * N));

for ii = 1:M
    ir = squeeze(IR(ii, :, :))';
    iririr = [ir; ir; ir];
    iririr2 = resample(iririr, targetFs, currentFs);
    N2 = round(length(iririr2)/3);
    ir2 = iririr2(N2+1:2*N2, :);
    IR2(ii, :, :) = ir2';
end

Obj.Data.IR = IR2;
Obj.Data.SamplingRate = targetFs;
Obj=SOFAupdateDimensions(Obj);

end


% Get aduio filenames routine
function audioFilenames = getAudioFilenames(params)

audioFilenames = cell(length(params.RecordingsScenarios), 1);

for iScenario = 1:length(params.RecordingsScenarios)
    scenario = params.RecordingsScenarios(iScenario);
    audioFilenames{iScenario} = dirWithoutDots(...
        fullfile(params.RecordingsBaseDir, scenario));
end

audioFilenames = cell2mat(audioFilenames);

end


% load tracks routine
function tracks = loadAudioTracks(audioFilename, params)

songName = fullfile(audioFilename.folder, audioFilename.name);
trackFilenames = dir(fullfile(songName, '*.wav'));
audioInfo = audioinfo(fullfile(trackFilenames(1).folder, ...
    trackFilenames(1).name));
totalSamples = audioInfo.TotalSamples ...
    - params.RecordingLoadRange(1) * params.RecordingsExpectedFs;
tracks = zeros(totalSamples, length(trackFilenames));

for iTrackFilename = 1:length(trackFilenames)
    trackPath = fullfile(trackFilenames(iTrackFilename).folder, ...
        trackFilenames(iTrackFilename).name);
    [track,Fs] = audioread(trackPath, ...
        params.RecordingLoadRange * params.RecordingsExpectedFs + [1 0]);
    
    if Fs ~= params.RecordingsExpectedFs
        error('Track frequency is not expected frequency');
    end
    
    tracks(:, iTrackFilename) = track;
end

end


% Spatialize all audio trakcs routine
% spatResults shape (width, HRTF, dir, sample, ch)
function [spatResults,spatMetaresults] = spatializeSong(HRTFs, tracks, params)

sz = [size(params.EnsambleAngleWidths, 2), ...
    length(HRTFs), ...
    size(params.EnsambleDirections, 1)];
dur = params.RecordingSpatRange(2) * params.RecordingsExpectedFs;
spatResults = zeros([sz dur params.NChannels]);
spatMetaresults = cell(sz);

for comb = allcomb(...
        1:length(params.EnsambleAngleWidths), ...
        1:length(HRTFs))'
    cComb = num2cell(comb);
    [iEnsambleAngleWidth,iHRTF] = cComb{:};
    
    ensambleAngleWidth = params.EnsambleAngleWidths(iEnsambleAngleWidth);
    
    for iEnsambleDirection = 1:size(params.EnsambleDirections, 1)
        ensambleDirection = params.EnsambleDirections(iEnsambleDirection, :);
        
        [spatResult,spatMetaresult] = spatializeAudioTracks(...
            tracks, HRTFs(iHRTF), ...
            ensambleDirection, ...
            ensambleAngleWidth, params);
        spatResults(iEnsambleAngleWidth,iHRTF,...
            iEnsambleDirection, :, :) = spatResult;
        spatMetaresults{iEnsambleAngleWidth,iHRTF,...
            iEnsambleDirection} = spatMetaresult;
        
        fprintf('Progress  [width %d/%d, HRTF %d/%d, dir %d/%d]\n', ...
            iEnsambleAngleWidth, ...
            size(params.EnsambleAngleWidths, 2), ...
            iHRTF, ...
            length(HRTFs), ...
            iEnsambleDirection, ...
            size(params.EnsambleDirections, 1));
    end
    
end

spatMetaresults = cell2mat(spatMetaresults);
end


% Spatialize audio routine
function [spatScene,spatMetaresult] = spatializeAudioTracks(...
    tracks, HRTF, ensambleDirection, ensambleAngleWidth, params)

spatMetaresult.RandTrackAngles = randSphCap(size(tracks, 2), ...
    ensambleDirection, ensambleAngleWidth);
desiredAnglePairs = [HRTF.Position];

[spatMetaresult.FittedHRTFAngles,spatMetaresult.FittedHRTFAnglesI] ...
    = findBestFitAnglePairs(spatMetaresult.RandTrackAngles, ...
    desiredAnglePairs, ensambleDirection, ensambleAngleWidth, ...
    params.EnsambleAngleWidthEps);

spatScene = [];

for iTrack = 1:size(tracks, 2)
    track = tracks(:, iTrack);
    iAngles = spatMetaresult.FittedHRTFAnglesI(iTrack);
    
    spatTrack = [
        conv(squeeze(HRTF.SOFA.Data.IR(iAngles, 1, :)), track) ...
        conv(squeeze(HRTF.SOFA.Data.IR(iAngles, 2, :)), track)];
    
%     azel = spatMetaresult.FittedHRTFAngles(iTrack, :);
%     spatTrack = SOFAspat(track, HRTF.SOFA, azel(1), azel(2));
    
    if isempty(spatScene)
        spatScene = zeros(size(spatTrack));
    end
    
    spatScene = spatScene + spatTrack;
end

spatScene = trimAndFadeSignal(spatScene, params);

end


% Trim and fade signal routine
function y = trimAndFadeSignal(x, params)

range = params.RecordingSpatRange * params.RecordingsExpectedFs - [0 1];
y = x(range(1):sum(range), :);

env = envGen(params.RecordingFadeTime(1), ...
    params.RecordingSpatRange(2), ...
    params.RecordingFadeTime(2), ...
    params.RecordingsExpectedFs, 2, 'sinsq')';
y = y .* env;

end


% Postprocess spatialization results routine
function spatResults = posprocessSpatResults(spatResults, params)

% Peak normalization and scaling
peakLevel = max(abs(spatResults), [], [3 4 5]);
spatResults = params.RecordingLevelScale * spatResults ./ peakLevel;

% DC equalization
spatResults = spatResults - mean(spatResults, 4);

end


% Save spatialization results routine
% spatResults shape (width, HRTF, dir, sample, ch)
function spatResults = saveSpatResults(spatResults, spatMetaresults, ...
    HRTFs, audioFilename, params)

if ~exist(params.SpatOutputDir, 'dir')
    mkdir(params.SpatOutputDir);
end

if ~exist(params.MetaOutputDir, 'dir')
    mkdir(params.MetaOutputDir);
end

for comb = allcomb(...
        1:length(params.EnsambleAngleWidths), ...
        1:length(HRTFs), ...
        1:size(params.EnsambleDirections, 1))'
    cComb = num2cell(comb);
    [iEnsambleAngleWidth,iHRTF,iEnsambleDirection] = cComb{:};
    
    spatFilename = getOutputFilename(...
        iEnsambleAngleWidth, iEnsambleDirection, HRTFs, iHRTF, ...
        audioFilename, params);
    
    scenario = getScenario(audioFilename);
    
    spatParentDir = fullfile(params.SpatOutputDir, ...
        num2str(params.EnsambleAngleWidths(iEnsambleAngleWidth)), ...
        scenario);
    metaParentDir = fullfile(params.MetaOutputDir, ...
        num2str(params.EnsambleAngleWidths(iEnsambleAngleWidth)), ...
        scenario);
    fullSpatFilename = fullfile(spatParentDir, spatFilename + '.wav');
    fullMetaFilename = fullfile(metaParentDir, ...
        spatFilename + '-spatMetaresults');
    
    if ~exist(spatParentDir, 'dir')
        mkdir(spatParentDir);
    end
    
    if ~exist(metaParentDir, 'dir')
        mkdir(metaParentDir);
    end
    
    spatOut = squeeze(spatResults(iEnsambleAngleWidth, iHRTF, ...
        iEnsambleDirection, :, :));
    audiowrite(fullSpatFilename, spatOut, ...
        params.RecordingsExpectedFs, ...
        'BitsPerSample', 32);
    save(fullMetaFilename, 'spatMetaresults', 'params', '-v7.3');
end

end


% Get scenario routine
function scenario = getScenario(audioFilename)

fileParts = split(audioFilename.folder, filesep);
scenario = fileParts{end};

end


% Get output filename routine
function [filename,parentDir] = getOutputFilename(iEnsambleAngleWidth, ...
    iEnsambleDirection, ...
    HRTFs, iHRTF, ...
    audioFilename, params)

HRTFGroup = HRTFs(iHRTF).HRTFGroup;
sceneName = params.EnsambleSceneNames(iEnsambleDirection);
ensambleDirection = params.EnsambleDirections(iEnsambleDirection, :);
ensambleAngleWidth = params.EnsambleAngleWidths(iEnsambleAngleWidth);

filename = sprintf("%s_hrtf%d_%s_scene%d_%s_az%d_el%d_width%d", ...
    audioFilename.name, ...
    iHRTF, ...
    HRTFGroup, ...
    iEnsambleDirection, ...
    sceneName, ...
    ensambleDirection(1), ...
    ensambleDirection(2), ...
    ensambleAngleWidth);

end


% Drawning routine
% spatMetaresults shape (width, HRTF, dir, audio)
function plotAudioScene(HRTFs, spatMetaresults, params)

if ~exist(params.FinalResultsOutputDir, 'dir')
    mkdir(params.FinalResultsOutputDir);
end

HRTFGroups = convertCharsToStrings(unique({HRTFs.HRTFGroup}));

for iHRTFGroup = 1:length(HRTFGroups)
    fig = figure('Position', [0, 0, 1400, 1200]);
    
    HRTFGroupName = HRTFGroups(iHRTFGroup);    
    HRTFidx = strcmp({HRTFs.HRTFGroup}, HRTFGroupName);
    
    for iEnsambleAngleWidth = 1:length(params.EnsambleAngleWidths)
        for iEnsambleDirection = 1:size(params.EnsambleDirections, 1)
            ensambleAngleWidth = params.EnsambleAngleWidths(iEnsambleAngleWidth);
            ensambleDirection = params.EnsambleDirections(iEnsambleDirection, :);

            selectedSpatMetaresults = spatMetaresults(...
                iEnsambleAngleWidth, HRTFidx, iEnsambleDirection, :);
            spatMetaresult = reshape(selectedSpatMetaresults, 1, []);

            fittedHRTFAngles = cat(1,spatMetaresult.FittedHRTFAngles);
            randTrackAngles = cat(1,spatMetaresult.RandTrackAngles);

            m = length(params.EnsambleAngleWidths);
            n = size(params.EnsambleDirections, 1);
            subplot(m, n, n * (iEnsambleAngleWidth - 1) + iEnsambleDirection);

            HRTFpos = unique(cat(1, HRTFs(HRTFidx).Position), 'rows');
            plotGroupedAnglePairs(HRTFpos, fittedHRTFAngles);

            [x, y, z] = sph2cart(...
                deg2rad(randTrackAngles(:, 1)), ...
                deg2rad(randTrackAngles(:, 2)), 1);

            hold on
            h = scatter3(x, y, z, 'filled');
            set(h, 'MarkerEdgeAlpha', 0.05, 'MarkerFaceAlpha', 0.05);
            hold off

            title(sprintf("width %d°, scene '%s' (%d°, %d°)", ...
                ensambleAngleWidth, ...
                params.EnsambleSceneNames(iEnsambleDirection), ...
                ensambleDirection(1), ensambleDirection(2)));
        end
    end
    
    sgtitle(HRTFGroupName);
    
    saveas(fig, fullfile(params.FinalResultsOutputDir, ...
        HRTFGroupName + '.png'));
end

end

% Comments
% params.EnsambleAngleWidths = [15 30];
% params.EnsambleDirections = [
%     0 0
%     -180 0];
% Elapsed time is 2471.161771 seconds.

% Progress  [audio: 141/152] (TalkToMeBaby)
% 23 HRTfs
% Elapsed time is 8941.610715 seconds.
% Elapsed time is 7079.953476 seconds. aftert changes in findBestFitAngle
% Elapsed time is 5176.862806 seconds. th-koln only
