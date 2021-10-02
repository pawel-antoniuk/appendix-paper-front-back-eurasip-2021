% Paweł Antoniuk 2021
% Bialystok University of Technology

function [timeArray, freqArray, specgram] = createSpecgram(inputSignal, params)
	maxDuration = params.AutioTotalSamples / params.AudioSampleRate;
    inputSignal = inputSignal(1:params.AudioDuration * params.AutioTotalSamples / maxDuration);
    
    key = params.SpecgramKey;
    
%     Create a spectrum array (spectrogram)
    powerHzScaleFactor = 0.5 * params.AudioSampleRate * sum(params.Window.^2);
    frames = v_enframe(inputSignal, params.Window, params.WindowHopSamples);
    spectrumArray = abs(v_rfft(frames, params.WindowLengthSamples, 2)) .^ 2 / powerHzScaleFactor;
    
%     Time params
    pSampleFrequency = params.AudioSampleRate / params.WindowHopSamples;
    pFirstSampleTime = 0.5 * (params.WindowLengthSamples + 1) / params.AudioSampleRate;
    pHop = params.AudioSampleRate / params.WindowLengthSamples;
    pFs = [pSampleFrequency, pFirstSampleTime, pHop];
    
%     Freqency params
    pFStep = (params.FHigh - params.FLow)/(params.NChannels - 1);
    pFRange = [params.FLow pFStep params.FHigh];
    
%     Create the final spectrogram using filter bank
    [timeArray, freqArray, specgram] = v_spgrambw(spectrumArray, pFs, ...
        key, 200, pFRange, params.DbRange);
end