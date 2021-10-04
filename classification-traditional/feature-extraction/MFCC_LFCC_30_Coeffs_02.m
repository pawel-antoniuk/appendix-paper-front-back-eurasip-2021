% Sławomir Zieliński 2021
% Bialystok University of Technology

%% Initialization
clear, close all, clc

experimentName = 'MFCC_LFCC_30_Coeffs';

inputFolderName = 'D:\FB Data\Binaural Audio Recordings';
outputFolderName = 'D:\FB Data\CSV';
outFileName = [outputFolderName,'\',experimentName,'_Features.csv']

Fs = 48000;

fdur = 0.02; % Frame duration in sec.
hop = 0.5;  % Hop factor
nfdur = fdur*Fs;
nhop = round(hop*nfdur);

f_low = 100;        % Low-frequency limit in Hz
f_high = 16000;     % High-frequency limit in Hz

nCoeff = 30;    % Number of cepstral coefficients (ncc=1 is the 0th coeff.)
nBands = 42;    % Number of channels in a filterbank

disp(['Number of MFCC features: ',num2str(2*4*nCoeff)]);
disp(['Number of LFCC features: ',num2str(2*4*nCoeff)]);


%% Filterbank
% Mel-frequency filterbank
melFreqFilterBank = designAuditoryFilterBank(Fs,...
    'FrequencyScale','mel',...
    'FFTLength',nfdur,...
    'FrequencyRange',[f_low,f_high],...
    'NumBands',nBands,...
    'Normalization','none',...
    'OneSided',true,...
    'FilterBankDesignDomain','linear');

% Custom lin-freq filterbank
numBandEdges = nBands + 2;
NFFT = nfdur;
filterBank = zeros(nBands,NFFT);
% bandEdges = logspace(log10(f_low),log10(f_high),numBandEdges);
bandEdges = linspace(f_low,f_high,numBandEdges);
bandEdgesBins = round((bandEdges/Fs)*NFFT) + 1;
for ii = 1:nBands
     filt = triang(bandEdgesBins(ii+2)-bandEdgesBins(ii));
     leftPad = bandEdgesBins(ii);
     rightPad = NFFT - numel(filt) - leftPad;
     filterBank(ii,:) = [zeros(1,leftPad),filt',zeros(1,rightPad)];
end
linFreqFilterBank = filterBank(:,1:NFFT/2+1);

% Plots
f_vector = 0.001*(0:floor(nfdur/2))*Fs/nfdur;
subplot(2,1,1)
plot(f_vector,melFreqFilterBank)
xlabel('Frequency (kHz)')
title('Mel-Frequency Filterbank Characteristics')
subplot(2,1,2)
plot(f_vector,linFreqFilterBank)
xlabel('Frequency (kHz)')
title('Lin-Frequency Filterbank Characteristics')



tStart = tic;


%% ======================= Feature Extraction =============================
disp('====== Feature Extraction ======')

%% Identify File Names
[fileNames, N] = getFileNames(inputFolderName);

%% Main Loop ==============================================================
Features = zeros(length(fileNames),2*2*4*nCoeff);
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
    
    %% Spectrum
    XX = stft(x,'Window',hamming(nfdur),'FFTLength',nfdur,'OverlapLength',nhop,'FrequencyRange','onesided');
    YY = stft(y,'Window',hamming(nfdur),'FFTLength',nfdur,'OverlapLength',nhop,'FrequencyRange','onesided');
    MM = stft(m,'Window',hamming(nfdur),'FFTLength',nfdur,'OverlapLength',nhop,'FrequencyRange','onesided');
    SS = stft(s,'Window',hamming(nfdur),'FFTLength',nfdur,'OverlapLength',nhop,'FrequencyRange','onesided');
 
    XX = abs(XX);
    YY = abs(YY);
    MM = abs(MM);
    SS = abs(SS);
    
    %% Cepstral coefficients
    x_MFCCs = cepstralCoefficients(melFreqFilterBank*XX,'NumCoeffs',nCoeff,'Rectification','log');
    y_MFCCs = cepstralCoefficients(melFreqFilterBank*YY,'NumCoeffs',nCoeff,'Rectification','log');
    m_MFCCs = cepstralCoefficients(melFreqFilterBank*MM,'NumCoeffs',nCoeff,'Rectification','log');
    s_MFCCs = cepstralCoefficients(melFreqFilterBank*SS,'NumCoeffs',nCoeff,'Rectification','log');

    x_LFCCs = cepstralCoefficients(linFreqFilterBank*XX,'NumCoeffs',nCoeff,'Rectification','log');
    y_LFCCs = cepstralCoefficients(linFreqFilterBank*YY,'NumCoeffs',nCoeff,'Rectification','log');
    m_LFCCs = cepstralCoefficients(linFreqFilterBank*MM,'NumCoeffs',nCoeff,'Rectification','log');
    s_LFCCs = cepstralCoefficients(linFreqFilterBank*SS,'NumCoeffs',nCoeff,'Rectification','log');    
    
    
    %% Summary statistics
    x_MFCCs_m = mean(x_MFCCs);
    y_MFCCs_m = mean(y_MFCCs);
    m_MFCCs_m = mean(m_MFCCs);
    s_MFCCs_m = mean(s_MFCCs);
    
    x_MFCCs_sd = std(x_MFCCs);
    y_MFCCs_sd = std(y_MFCCs);
    m_MFCCs_sd = std(m_MFCCs);
    s_MFCCs_sd = std(s_MFCCs); 
    
    x_LFCCs_m = mean(x_LFCCs);
    y_LFCCs_m = mean(y_LFCCs);
    m_LFCCs_m = mean(m_LFCCs);
    s_LFCCs_m = mean(s_LFCCs);
    
    x_LFCCs_sd = std(x_LFCCs);
    y_LFCCs_sd = std(y_LFCCs);
    m_LFCCs_sd = std(m_LFCCs);
    s_LFCCs_sd = std(s_LFCCs);     
       
    Features(ii, :) = [ x_MFCCs_m ...
            y_MFCCs_m ...
            m_MFCCs_m ...
            s_MFCCs_m ...
            x_MFCCs_sd ...
            y_MFCCs_sd ...
            m_MFCCs_sd ...
            s_MFCCs_sd ...
            x_LFCCs_m ...
            y_LFCCs_m ...
            m_LFCCs_m ...
            s_LFCCs_m ...
            x_LFCCs_sd ...
            y_LFCCs_sd ...
            m_LFCCs_sd ...
            s_LFCCs_sd ...            
          ];      
    
end % End of loop across files 


%% Variable names
vars_MFCC_x_m = sequenceofVariables('mfcc_x_m_',nCoeff);
vars_MFCC_y_m = sequenceofVariables('mfcc_y_m_',nCoeff);
vars_MFCC_m_m = sequenceofVariables('mfcc_m_m_',nCoeff);
vars_MFCC_s_m = sequenceofVariables('mfcc_s_m_',nCoeff);

vars_MFCC_x_sd = sequenceofVariables('mfcc_x_sd_',nCoeff);
vars_MFCC_y_sd = sequenceofVariables('mfcc_y_sd_',nCoeff);
vars_MFCC_m_sd = sequenceofVariables('mfcc_m_sd_',nCoeff);
vars_MFCC_s_sd = sequenceofVariables('mfcc_s_sd_',nCoeff);

vars_LFCC_x_m = sequenceofVariables('lfcc_x_m_',nCoeff);
vars_LFCC_y_m = sequenceofVariables('lfcc_y_m_',nCoeff);
vars_LFCC_m_m = sequenceofVariables('lfcc_m_m_',nCoeff);
vars_LFCC_s_m = sequenceofVariables('lfcc_s_m_',nCoeff);

vars_LFCC_x_sd = sequenceofVariables('lfcc_x_sd_',nCoeff);
vars_LFCC_y_sd = sequenceofVariables('lfcc_y_sd_',nCoeff);
vars_LFCC_m_sd = sequenceofVariables('lfcc_m_sd_',nCoeff);
vars_LFCC_s_sd = sequenceofVariables('lfcc_s_sd_',nCoeff);

variableNames = [vars_MFCC_x_m vars_MFCC_y_m vars_MFCC_m_m vars_MFCC_s_m...
    vars_MFCC_x_sd vars_MFCC_y_sd vars_MFCC_m_sd vars_MFCC_s_sd...
    vars_LFCC_x_m vars_LFCC_y_m vars_LFCC_m_m vars_LFCC_s_m...
    vars_LFCC_x_sd vars_LFCC_y_sd vars_LFCC_m_sd vars_LFCC_s_sd];


% Create table
dt = array2table(Features, 'VariableNames',variableNames);
dt = addvars(dt,fileNames,'Before',1);


%% Save data
outputfile = outFileName;
writetable(dt,outputfile);


%% Finalize
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));











