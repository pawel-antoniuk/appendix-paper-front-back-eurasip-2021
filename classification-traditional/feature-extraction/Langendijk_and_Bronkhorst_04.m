%% Initialization
clear, close all, clc

experimentName = 'Langendijk_and_Bronkhorst';

inputFolderName = 'D:\FB Data\Binaural Audio Recordings';
outputFolderName = 'D:\FB Data\CSV';
outFileName = [outputFolderName,'\',experimentName,'_Features.csv'];

Fs = 48000;
filterOrder = 512;

B1_f1 = 8000;
B1_f2 = 16000;


fdur = 0.02; % Frame duration in sec.
hop = 0.5;  % Hop factor
nfdur = fdur*Fs;
nhop = round(hop*nfdur);

audioFE = audioFeatureExtractor( ...
'SampleRate',Fs, ...
'Window',hamming(round(nfdur),'periodic'), ...
'OverlapLength',round(nhop), ...
'SpectralDescriptorInput','linearSpectrum', ...
'spectralCentroid',true, ...
'spectralCrest',true, ...
'spectralDecrease',true, ...
'spectralEntropy',true, ...
'spectralFlatness',true, ...
'spectralFlux',true, ...
'spectralKurtosis',true, ...
'spectralRolloffPoint',true, ...
'spectralSkewness',true, ...
'spectralSlope',true, ...
'spectralSpread',true);

info(audioFE)


%% Filter design
filter1 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B1_f1,'CutoffFrequency2',B1_f2, ...
         'SampleRate',Fs);

     
% Plot filter characteristic
% fvtool(filter1);


tStart = tic;


%% ======================= Feature Extraction =============================
disp('====== Feature Extraction ======')


%% Identify File Names
[fileNames, N] = getFileNames(inputFolderName);


%% Main Loop ==============================================================
Features = zeros(length(fileNames), 11*4*2);
parfor ii = 1:N % Loop across files

    fname = fullfile(inputFolderName,fileNames{ii});    
%     disp(' ' );
%     disp(['Progress: ',num2str(ii),' out of ',num2str(N)])
%     disp(fname);

    %% Open Audio File
    [xy, fs] = audioread(fname);
    xy = xy - ones(size(xy))*diag(mean(xy)); % Remove DC offset  
    xy = xy / max(rms(xy)); % RMS joint-channel normalization
    
    x = xy(:,1);
    y = xy(:,2);

    %% Mid-Side Processing   
    m = xy * [1; 1]; % m-signal
    s = xy * [1; -1]; % s-signal
    m = m - mean(m); % Remove DC offset
    s = s - mean(s); % Remove DC offset

    %% Filtering
    x1 = filtfilt(filter1, x);   
    y1 = filtfilt(filter1, y); 
    m1 = filtfilt(filter1, m);
    s1 = filtfilt(filter1, s);
    
    
    %% Plot
    % subplot(2,1,1); pwelch(x,2048,1024,[1:20000],Fs);
    % subplot(2,1,2); pwelch(x1,2048,1024,[1:20000],Fs); 
    
    
    %% Audio Features   
    x1_featrures = extract(audioFE,x1);
    y1_featrures = extract(audioFE,y1);
    m1_featrures = extract(audioFE,m1);
    s1_featrures = extract(audioFE,s1);

    % Plot
    % plot(x1_featrures)
    
    
        %% Summary statistics
    % Mean values
    x1_m = mean(x1_featrures);   
    y1_m = mean(y1_featrures);
    m1_m = mean(m1_featrures);
    s1_m = mean(s1_featrures);

    
    % Standard deviations
    x1_sd = std(x1_featrures);
    y1_sd = std(y1_featrures);
    m1_sd = std(m1_featrures);
    s1_sd = std(s1_featrures);
 
    
    % Features
    Features(ii, :) = [ x1_m ...
                        y1_m ...
                        m1_m ...
                        s1_m ...
                        x1_sd ...
                        y1_sd ...
                        m1_sd ...
                        s1_sd ...                               
                        ];
                  
    
end % End of loop across files 


%% Listen
% clip = audioplayer(0.5*x1, Fs);      
% play(clip);   
% % stop(clip)


%% Variable names
variableNames = {
    'centroid_x_m'
    'crest_x_m'
    'decrease_x_m'
    'entropy_x_m'
    'flatness_x_m'
    'flux_x_m'
    'kurtosis_x_m'
    'rolloff_x_m'
    'skewness_x_m'
    'slope_x_m'
    'spread_x_m'
    
    'centroid_y_m'
    'crest_y_m'
    'decrease_y_m'
    'entropy_y_m'
    'flatness_y_m'
    'flux_y_m'
    'kurtosis_y_m'
    'rolloff_y_m'
    'skewness_y_m'
    'slope_y_m'
    'spread_y_m'
    
    'centroid_m_m'
    'crest_m_m'
    'decrease_m_m'
    'entropy_m_m'
    'flatness_m_m'
    'flux_m_m'
    'kurtosis_m_m'
    'rolloff_m_m'
    'skewness_m_m'
    'slope_m_m'
    'spread_m_m'
    
    'centroid_s_m'
    'crest_s_m'
    'decrease_s_m'
    'entropy_s_m'
    'flatness_s_m'
    'flux_s_m'
    'kurtosis_s_m'
    'rolloff_s_m'
    'skewness_s_m'
    'slope_s_m'
    'spread_s_m'
    
    'centroid_x_sd'
    'crest_x_sd'
    'decrease_x_sd'
    'entropy_x_sd'
    'flatness_x_sd'
    'flux_x_sd'
    'kurtosis_x_sd'
    'rolloff_x_sd'
    'skewness_x_sd'
    'slope_x_sd'
    'spread_x_sd'
    
    'centroid_y_sd'
    'crest_y_sd'
    'decrease_y_sd'
    'entropy_y_sd'
    'flatness_y_sd'
    'flux_y_sd'
    'kurtosis_y_sd'
    'rolloff_y_sd'
    'skewness_y_sd'
    'slope_y_sd'
    'spread_y_sd'
    
    'centroid_m_sd'
    'crest_m_sd'
    'decrease_m_sd'
    'entropy_m_sd'
    'flatness_m_sd'
    'flux_m_sd'
    'kurtosis_m_sd'
    'rolloff_m_sd'
    'skewness_m_sd'
    'slope_m_sd'
    'spread_m_sd'
    
    'centroid_s_sd'
    'crest_s_sd'
    'decrease_s_sd'
    'entropy_s_sd'
    'flatness_s_sd'
    'flux_s_sd'
    'kurtosis_s_sd'
    'rolloff_s_sd'
    'skewness_s_sd'
    'slope_s_sd'
    'spread_s_sd'    
    };


% Create table
dt = array2table(Features, 'VariableNames',variableNames);
dt = addvars(dt,fileNames,'Before',1);


%% Save data
outputfile = outFileName;
writetable(dt,outputfile);


%% Finalize
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));













