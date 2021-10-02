%% Initialization
clear, close all, clc

experimentName = 'Blauert_Boosted_Bands';

inputFolderName = 'D:\FB Data\Binaural Audio Recordings';
outputFolderName = 'D:\FB Data\CSV';
outFileName = [outputFolderName,'\',experimentName,'_Features.csv'];

Fs = 48000;
filterOrder = 512;

B1_f1 = 150;
B1_f2 = 540;

B2_f1 = 720;
B2_f2 = 1700;

B3_f1 = 1900;
B3_f2 = 2900;

B4_f1 = 3600;
B4_f2 = 5800;

B5_f1 = 7400;
B5_f2 = 11100;

fdur = 0.02; % Frame duration in sec.
hop = 0.5;  % Hop factor
nfdur = fdur*Fs;
nhop = round(hop*nfdur);



%% Filter design
filter1 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B1_f1,'CutoffFrequency2',B1_f2, ...
         'SampleRate',Fs);
filter2 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B2_f1,'CutoffFrequency2',B2_f2, ...
         'SampleRate',Fs);
filter3 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B3_f1,'CutoffFrequency2',B3_f2, ...
         'SampleRate',Fs);
filter4 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B4_f1,'CutoffFrequency2',B4_f2, ...
         'SampleRate',Fs);
filter5 = designfilt('bandpassfir','FilterOrder',filterOrder, ...
         'CutoffFrequency1',B5_f1,'CutoffFrequency2',B5_f2, ...
         'SampleRate',Fs);     
     
% Plot filter characteristic
% fvtool(filter5);



tStart = tic;


%% ======================= Feature Extraction =============================
disp('====== Feature Extraction ======')


%% Identify File Names
[fileNames, N] = getFileNames(inputFolderName);

%% Main Loop ==============================================================
Features = zeros(length(fileNames), 5*4*2);
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
    x2 = filtfilt(filter2, x);
    x3 = filtfilt(filter3, x);
    x4 = filtfilt(filter4, x);
    x5 = filtfilt(filter5, x);    
    
    y1 = filtfilt(filter1, y);
    y2 = filtfilt(filter2, y);
    y3 = filtfilt(filter3, y);
    y4 = filtfilt(filter4, y);
    y5 = filtfilt(filter5, y); 
    
    m1 = filtfilt(filter1, m);
    m2 = filtfilt(filter2, m);
    m3 = filtfilt(filter3, m);
    m4 = filtfilt(filter4, m);
    m5 = filtfilt(filter5, m);  
    
    s1 = filtfilt(filter1, s);
    s2 = filtfilt(filter2, s);
    s3 = filtfilt(filter3, s);
    s4 = filtfilt(filter4, s);
    s5 = filtfilt(filter5, s);       


    %% Plot
    % subplot(2,1,1); pwelch(s,2048,1024,[1:20000],Fs);
    % subplot(2,1,2); pwelch(s3,2048,1024,[1:20000],Fs);    
    
    %% Frame-based Processing using Voicebox    
    [X1,tfrm]=v_enframe(x1,hamming(nfdur),nhop);
    [X2,~]=v_enframe(x2,hamming(nfdur),nhop);
    [X3,~]=v_enframe(x3,hamming(nfdur),nhop);
    [X4,~]=v_enframe(x4,hamming(nfdur),nhop);
    [X5,~]=v_enframe(x5,hamming(nfdur),nhop);    
    
    [Y1,~]=v_enframe(y1,hamming(nfdur),nhop);
    [Y2,~]=v_enframe(y2,hamming(nfdur),nhop);
    [Y3,~]=v_enframe(y3,hamming(nfdur),nhop);
    [Y4,~]=v_enframe(y4,hamming(nfdur),nhop);
    [Y5,~]=v_enframe(y5,hamming(nfdur),nhop);  
    
    [M1,~]=v_enframe(m1,hamming(nfdur),nhop);
    [M2,~]=v_enframe(m2,hamming(nfdur),nhop);
    [M3,~]=v_enframe(m3,hamming(nfdur),nhop);
    [M4,~]=v_enframe(m4,hamming(nfdur),nhop);
    [M5,~]=v_enframe(m5,hamming(nfdur),nhop);     

    [S1,~]=v_enframe(s1,hamming(nfdur),nhop);
    [S2,~]=v_enframe(s2,hamming(nfdur),nhop);
    [S3,~]=v_enframe(s3,hamming(nfdur),nhop);
    [S4,~]=v_enframe(s4,hamming(nfdur),nhop);
    [S5,~]=v_enframe(s5,hamming(nfdur),nhop);    
    
    
    %% RMS
    x1_rms = 20*log10(rms(X1'));
    x2_rms = 20*log10(rms(X2'));
    x3_rms = 20*log10(rms(X3'));
    x4_rms = 20*log10(rms(X4'));
    x5_rms = 20*log10(rms(X5'));
    
    y1_rms = 20*log10(rms(Y1'));
    y2_rms = 20*log10(rms(Y2'));
    y3_rms = 20*log10(rms(Y3'));
    y4_rms = 20*log10(rms(Y4'));
    y5_rms = 20*log10(rms(Y5'));    
    
    m1_rms = 20*log10(rms(M1'));
    m2_rms = 20*log10(rms(M2'));
    m3_rms = 20*log10(rms(M3'));
    m4_rms = 20*log10(rms(M4'));
    m5_rms = 20*log10(rms(M5'));      
    
    s1_rms = 20*log10(rms(S1'));
    s2_rms = 20*log10(rms(S2'));
    s3_rms = 20*log10(rms(S3'));
    s4_rms = 20*log10(rms(S4'));
    s5_rms = 20*log10(rms(S5'));     
        
    % Plot
    % plot(tfrm/Fs,m5_rms)

    
    %% Summary statistics
    % Mean values
    x1_m = mean(x1_rms);
    x2_m = mean(x2_rms); 
    x3_m = mean(x3_rms);
    x4_m = mean(x4_rms);
    x5_m = mean(x5_rms);    
    
    y1_m = mean(y1_rms);
    y2_m = mean(y2_rms); 
    y3_m = mean(y3_rms);
    y4_m = mean(y4_rms);
    y5_m = mean(y5_rms);  
    
    m1_m = mean(m1_rms);
    m2_m = mean(m2_rms); 
    m3_m = mean(m3_rms);
    m4_m = mean(m4_rms);
    m5_m = mean(m5_rms);    
    
    s1_m = mean(s1_rms);
    s2_m = mean(s2_rms); 
    s3_m = mean(s3_rms);
    s4_m = mean(s4_rms);
    s5_m = mean(s5_rms);
    
    % Standard deviations
    x1_sd = std(x1_rms);
    x2_sd = std(x2_rms); 
    x3_sd = std(x3_rms);
    x4_sd = std(x4_rms);
    x5_sd = std(x5_rms);    
    
    y1_sd = std(y1_rms);
    y2_sd = std(y2_rms); 
    y3_sd = std(y3_rms);
    y4_sd = std(y4_rms);
    y5_sd = std(y5_rms);  
    
    m1_sd = std(m1_rms);
    m2_sd = std(m2_rms); 
    m3_sd = std(m3_rms);
    m4_sd = std(m4_rms);
    m5_sd = std(m5_rms);    
    
    s1_sd = std(s1_rms);
    s2_sd = std(s2_rms); 
    s3_sd = std(s3_rms);
    s4_sd = std(s4_rms);
    s5_sd = std(s5_rms);      
                   
    % Features
    Features(ii, :) = [ x1_m ...
                        x2_m ... 
                        x3_m ...
                        x4_m ...
                        x5_m ...   
                        y1_m ...
                        y2_m ... 
                        y3_m ...
                        y4_m ...
                        y5_m ...
                        m1_m ...
                        m2_m ...
                        m3_m ...
                        m4_m ...
                        m5_m ...    
                        s1_m ...
                        s2_m ...
                        s3_m ...
                        s4_m ...
                        s5_m ...
                        x1_sd ...
                        x2_sd ...
                        x3_sd ...
                        x4_sd ...
                        x5_sd ...   
                        y1_sd ...
                        y2_sd ...
                        y3_sd ...
                        y4_sd ...
                        y5_sd ... 
                        m1_sd ...
                        m2_sd ...
                        m3_sd ...
                        m4_sd ...
                        m5_sd ...  
                        s1_sd ...
                        s2_sd ... 
                        s3_sd ...
                        s4_sd ...
                        s5_sd ...                                                          
                        ];
    
end % End of loop across files 


%% Listen
% clip = audioplayer(0.5*s5, Fs);      
% play(clip);   
% % stop(clip


%% Variable names
vars_bb_x_m = sequenceofVariables('bb_x_m_',5);
vars_bb_y_m = sequenceofVariables('bb_y_m_',5);
vars_bb_mid_m = sequenceofVariables('bb_mid_m_',5);
vars_bb_side_m = sequenceofVariables('bb_side_m_',5);

vars_bb_x_sd = sequenceofVariables('bb_x_sd_',5);
vars_bb_y_sd = sequenceofVariables('bb_y_sd_',5);
vars_bb_mid_sd = sequenceofVariables('bb_mid_sd_',5);
vars_bb_side_sd = sequenceofVariables('bb_side_sd_',5);

variableNames = [vars_bb_x_m...
    vars_bb_y_m...
    vars_bb_mid_m...
    vars_bb_side_m...
    vars_bb_x_sd...
    vars_bb_y_sd...
    vars_bb_mid_sd...
    vars_bb_side_sd...
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












