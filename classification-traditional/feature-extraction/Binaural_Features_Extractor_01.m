%% Initialization
startTwoEars
clear, close all, clc

experimentName = 'Binaural';

inputFolderName = 'D:\FB Data\Binaural Audio Recordings';
outputFolderName = 'D:\FB Data\CSV';
outFileName = [outputFolderName,'\',experimentName,'_Features.csv'];

Fs = 48000;

fdur = 0.02; % Frame duration in sec.
hop = 0.5;  % Hop factor
nfdur = fdur*Fs;
nhop = round(hop*nfdur);

nch = 42;   % Number of channels in a filterbank



tStart = tic;

%% ======================= Feature Extraction =============================
disp('====== Feature Extraction ======')


%% Identify File Names
[fileNames, N] = getFileNames(inputFolderName);


%% Main Loop ==============================================================
Features = zeros(length(fileNames), nch*3*2);
parfor ii = 1:N % Loop across files

    fname = fullfile(inputFolderName,fileNames{ii});
%     disp(' ' );
%     disp(['Progress: ',num2str(ii),' out of ',num2str(N)])
%     disp(fname);

    %% Open Audio File
    [xy, fs] = audioread(fname);
    xy = xy - ones(size(xy))*diag(mean(xy)); % Remove DC offset  
    xy = xy / max(rms(xy)); % RMS joint-channel normalization
    
    % x = xy(:,1);
    % y = xy(:,2);
    
    dur = length(xy)/Fs; % Duration of a signal to be analyzed

    %% Mid-Side Processing   
    m = xy * [1; 1]; % m-signal
    s = xy * [1; -1]; % s-signal
    m = m - mean(m); % Remove DC offset
    s = s - mean(s); % Remove DC offset
    ms = [m s];

    %% Binaural Cues - Two!Ears
    
    % Parameters
    parameters = genParStruct(...
        'pp_bRemoveDC',true,...
        'pp_cutoffHzDC',20,...
        'fb_type','gammatone',...
        'fb_nChannels',nch,...
        'fb_lowFreqHz',100,...
        'fb_highFreqHz',16000,...
        'ihc_method','dau',...
        'ild_wSizeSec',fdur,...
        'ild_hSizeSec',hop*fdur,...
        'ild_wname','hann',...
        'cc_wSizeSec',fdur,...
        'cc_hSizeSec',hop*fdur,...
        'cc_wname','hann',...
        'cc_maxDelaySec',0.0011,...
        'rm_wSizeSec',fdur,...
        'rm_hSizeSec',hop*fdur,...
        'rm_scaling','power',...
        'rm_decaySec',8E-3,...
        'rm_wname','hann'); 

        % Requests
        requests = {'ild','itd','ic','spectralFeatures'}; 

        % Instantiation of a data and manager objects
        dObj = dataObject(xy,Fs,dur,2);
        dObj_ms = dataObject(ms,Fs,dur,2);
        mObj = manager(dObj, requests, parameters);
        mObj_ms = manager(dObj_ms, requests, parameters);

        % Request processing
        mObj.processSignal();
        mObj_ms.processSignal(); 
        
%         % Find out dependencies
%         Processor.getDependencyList('ildProc')
%         Processor.getDependencyList('itdProc')
%         parameterHelper % Use it to check manually for parameters          

        % Extracting the data
        ILDs = dObj.ild{1,1}.Data(:)';
        ITDs = dObj.itd{1,1}.Data(:)';
        ICs = dObj.ic{1,1}.Data(:)';

%         % Plot
%         figure
%         imagesc(ITDs); set(gca,'Ydir','Normal'); colorbar; title('Feature')

%         % Plots
%         sf_no = 14
%         subplot(2,4,1); plot(SF_x(sf_no,:))
%         subplot(2,4,2); plot(SF_y(sf_no,:))
%         subplot(2,4,3); plot(SF_m(sf_no,:))
%         subplot(2,4,4); plot(SF_s(sf_no,:))
%         subplot(2,4,5); plot(d_SF_x(sf_no,:))
%         subplot(2,4,6); plot(d_SF_y(sf_no,:))
%         subplot(2,4,7); plot(d_SF_m(sf_no,:))
%         subplot(2,4,8); plot(d_SF_s(sf_no,:))


        %% Summary statistics
        ILDs_m = mean(ILDs,2)';
        ILDs_sd = std(ILDs,0,2)';
        ITDs_m = mean(ITDs,2)';
        ITDs_sd = std(ITDs,0,2)';
        ICs_m = mean(ICs,2)';
        ICs_sd = std(ICs,0,2)';

        
        Features(ii, :) = [ ILDs_m ...
                            ILDs_sd ...
                            ITDs_m ...
                            ITDs_sd ...
                            ICs_m ...
                            ICs_sd ...
                            ];  

end % End of loop across files 

%% Variable names
vars_ild_m = sequenceofVariables('ild_m_ch_',nch);
vars_ild_sd = sequenceofVariables('ild_sd_ch_',nch);
vars_itd_m = sequenceofVariables('itd_m_ch_',nch);
vars_itd_sd = sequenceofVariables('itd_sd_ch_',nch);
vars_ic_m = sequenceofVariables('ic_m_ch_',nch);
vars_ic_sd = sequenceofVariables('ic_sd_ch_',nch);

variableNames = [vars_ild_m...
    vars_ild_sd...
    vars_itd_m...
    vars_itd_sd...
    vars_ic_m...
    vars_ic_sd...
    ];

% Create table
dt = array2table(Features, 'VariableNames',variableNames);
dt = addvars(dt,fileNames,'Before',1);


%% Save data
outputfile = outFileName;
writetable(dt,outputfile);

%% Finalize
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

