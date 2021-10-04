% Sławomir Zieliński 2021
% Bialystok University of Technology

%% HRTF-Based Binaural Audio Synthesis

% SOFA package required. 
SOFAstart
SOFAgetVersion()

%% Initialization
clear all, close all, clc

typeOfStimuli = 'Train';    % 'Train' or 'Test'

hrtfDir = 'D:\FB Data\HRTF Database';
inputAudioFolder = 'D:\FB Data\Multitrack Recordings\RMSNorm 48kHz';
outputAudioFolder = 'D:\FB Data\Binaural Audio Recordings';

NofScenes = 2;      % Front, Back
augmentationNumber = 1

viewingAngle = 60;  % Viewing angle for front scene, default value = 60 (180 for full hemifield)
targetFs = 48000;   % Target sample rate in Hz
ele = 0;            % Elevation
loadGap = 2.5;      % Time in sec. to skip whan audio loading
initGap = 0.5;      % Time in sec. to skip at the beginning during trimming
dur = 7;            % Target duration in sec. during trimming
fadeTime = 0.010;   % Fade-in and fade-out time in sec. during trimming
audioScaling = 0.9; % Scaling of saved audio files


%% Loading HRTF Datasets
disp('Loading HRTFs...');
hrtfNames = hrtfFileNames;

HRTF = struct([]);
for hh = 1:length(hrtfNames)
    disp(hrtfNames{hh});
    HRTF{hh} = SOFAload(fullfile(hrtfDir,hrtfNames{hh}));
    if HRTF{hh}.Data.SamplingRate ~= 48000
        error('Error! Sample rate different than 48kHz.')
    end
    % If the number of samples is odd, remove the last sample
    % It fixes the problem with SOFA-based convolution
    if mod(HRTF{hh}.API.N,2) ~= 0
        tmpIR = HRTF{hh}.Data.IR(:,:,1:end-1); % Remove last sample
        HRTF{hh}.Data.IR = tmpIR;
        HRTF{hh}.API.N = size(tmpIR,3);
    end        
end


%% Identify Folders with Multitrack Audio Recordings
inputAudioFolder = fullfile(inputAudioFolder,typeOfStimuli);
[folderNames, M] = getFolderNames(inputAudioFolder);
disp(['Number of folders = ', num2str(M)]);


%% Main Loops =============================================================
tStart = tic;
for hrtf_no = 1:length(hrtfNames) % Loop across HRTFs

    % Verify HRTF
    Obj = HRTF{hrtf_no}
    % SOFAinfo(Obj);
    size(Obj.Data.IR)
    disp(['min azimuth:',num2str(min(Obj.SourcePosition(:,1)))])
    disp(['max azimuth:',num2str(max(Obj.SourcePosition(:,1)))])
    disp(['min elevation:',num2str(min(Obj.SourcePosition(:,2)))])
    disp(['max elevation:',num2str(max(Obj.SourcePosition(:,2)))])

    parfor mm = 1:M % Loop across folders (multitrack recordings)
        disp(' ');
        disp(['HRTF ',num2str(hrtf_no)]);
        disp(['Progress: ',num2str(mm),' out of ',num2str(M)])
        
        %% Identify file names in current folder
        currentFolderName = fullfile(inputAudioFolder,folderNames{mm});
        disp(['Current folder:  ',currentFolderName]);
        [fileNames, N] = getFileNames(currentFolderName);
        disp(['Convolving ', num2str(N), ' files...']);      
        
        %% Loading Audio Files
        X = [];
        for ii = 1:N
            sourceFname = fullfile(currentFolderName,fileNames{ii});
            [x, fs] = audioread(sourceFname,[loadGap*targetFs+1,inf]);
            X(:,ii) = x;
            if fs ~= targetFs
                error('Sampling rate mismatch.');
            end            
        end  
    
        
        %% Randomly Choose Azimuth Angles
        aziFront = randomAngleSelector(N,-viewingAngle/2,viewingAngle/2); % Front sector
        
        % For back sector symmetricly flip front sources 
        % aziBack = -1 * aziFront + 180; 
        aziBack = zeros(1,length(aziFront));
        for aa = 1:length(aziFront)
            if aziFront(aa) >= 0 
                aziBack(aa) = 180 - aziFront(aa);
            else 
                aziBack(aa) = -180 - aziFront(aa);
            end
        end


        %% Plot angles
        close all
        txtTitle = [folderNames{mm},'      N = ',num2str(N)];
        plotAngles(aziFront,aziBack,txtTitle);
        drawnow 
        

        %% Convolution based on SOFAspat
        individuallyPannedSourcesFront = [];
        individuallyPannedSourcesBack = []; 
        for nn = 1:N
            individuallyPannedSourcesFront(nn,:,:) = ...
                SOFAspat(X(:,nn),Obj,aziFront(nn),ele);
            individuallyPannedSourcesBack(nn,:,:) = ...
                SOFAspat(X(:,nn),Obj,aziBack(nn),ele);
        end
        % Create scenes - Mixing individual signals (adding up)
        yFront = squeeze(sum(individuallyPannedSourcesFront,1));
        yBack = squeeze(sum(individuallyPannedSourcesBack,1));
               
      
        
        %% Control listening
        
%         % Listen to individual sources from the front
%         sourceID = 2
%         y = squeeze(individuallyPannedSourcesFront(sourceID,:,:));
%         y = 0.05*y/rms(y(:));
%         clip = audioplayer(y, targetFs);      
%         play(clip);   
%         % stop(clip)
%         disp(['Azimuth: ', num2str(aziFront(sourceID))])
%     
%         
%         % Listen to individual sources from the back
%         y = squeeze(individuallyPannedSourcesBack(sourceID,:,:));
%         y = 0.05*y/rms(y(:));
%         clip = audioplayer(y, targetFs);      
%         play(clip);   
%         % stop(clip)
%         disp(['Azimuth: ', num2str(aziBack(sourceID))])
%     
%         
%         % Listen to a composite audio scene - Front
%         y = yFront;
%         y = 0.05*y/rms(y(:));
%         max(abs(y(:)))
%         clip = audioplayer(y, targetFs);      
%         play(clip); 
%         % stop(clip)   
%         
%         
%         % Listen to a composite audio scene - Back
%         y = yBack;
%         y = 0.05*y/rms(y(:));
%         max(abs(y(:)))
%         clip = audioplayer(y, targetFs);      
%         play(clip); 
%         % stop(clip)        

        
        

        %% Trim with fade in and fade out
        yFront2 = yFront(initGap*targetFs:(initGap+dur)*targetFs-1,:);
        yBack2 = yBack(initGap*targetFs:(initGap+dur)*targetFs-1,:);
        env = envGen(fadeTime,dur,fadeTime,targetFs,2,"sinsq")';
        yFront3 = yFront2.*env;
        yBack3 = yBack2.*env;


        %% RMS normalization, scaling and DC equalization
%         % RMS joint-normalization
%         yFront4 = audioScaling*yFront3/rms(yFront3(:)); 
%         yBack4 = audioScaling*yBack3/rms(yBack3(:)); 
%         % DC equalization
%         yFront5 = yFront4 - ones(size(yFront4))*diag(mean(yFront4));
%         yBack5 = yBack4 - ones(size(yBack4))*diag(mean(yBack4));
%         % Safeguard against clipping
%         peakLevel = max(abs([yFront4(:); yBack4(:)]));
%         if peakLevel >= 1
%             disp(['PeakLevel = ',num2str(peakLevel)]);
%             error('Signal(s) clipped!');
%         end
%         Y = zeros(NofScenes,size(yFront5,1),2);
%         Y(1,:,:) = yFront5;
%         Y(2,:,:) = yBack5;    
        

        %% Peak normalization, scaling and DC equalization
        peakLevel = max(abs([yFront3(:); yBack3(:)]));
        yFront4 = audioScaling*yFront3/peakLevel; 
        yBack4 = audioScaling*yBack3/peakLevel; 
        % DC equalization
        yFront5 = yFront4 - ones(size(yFront4))*diag(mean(yFront4));
        yBack5 = yBack4 - ones(size(yBack4))*diag(mean(yBack4)); 
        Y = zeros(NofScenes,size(yFront5,1),2);
        Y(1,:,:) = yFront5;
        Y(2,:,:) = yBack5;


        %% Saving
        hrtf_name_split = split(hrtfNames{hrtf_no},'\');
        hrtf_source = hrtf_name_split{1};
        for sc = 1:NofScenes
            fileName = [folderNames{mm},'_hrtf',num2str(hrtf_no),...
                '_',hrtf_source, ...
                '_scene',num2str(sc),'_',sceneDecoder(sc),...
                '_aug',num2str(augmentationNumber),'.wav'];      
            fullFname = fullfile(outputAudioFolder,typeOfStimuli,fileName);
            audiowrite(fullFname,squeeze(Y(sc,:,:)),targetFs,'BitsPerSample',32);
        end         
        
        
    end % End of loop across folders

end % End of loop across HRTFs

%% Finalize
disp(' ')
disp('Convolution finished.')
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));



%% Explore selected HRTF
% Obj = HRTF{hrtf_no}
% SOFAdetails(Obj)
% size(Obj.Data.IR)
% SOFAplotGeometry_SZ(Obj)
% ASV = SOFAcalculateAPV(Obj); % Apparent Source Vector
% elev_res = unique(sort(ASV(:,2))) % elevation angles
% azim_res = unique(sort(ASV(:,1))) % azim angles
% % SOFAinfo(Obj)






