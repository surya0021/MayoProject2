% This program gets measures used for single electrode analysis. The
% following measures are used:
% 1. Spike Counts
% 2. LFP Power using MT
% 3. LFP amplitude and phase using FFT

function [allFiringRates,allMTPower,allFFTVals,allNumTrials,allTargetOnsetTimes,freqValsMT,freqValsFFT] = getAnalysisMeasuresSingleElectrode(TWNum)

if ~exist('TWNum','var');                   TWNum=3;                    end

tapers = [TWNum 2*TWNum-1];

folderSavedData = fullfile(pwd,'savedData');
fileNameSave = fullfile(folderSavedData,['singleElectrodeMeasures' num2str(TWNum) '.mat']);

if exist(fileNameSave,'file')
    disp(['Loading saved data in ' fileNameSave]);
    load(fileNameSave,'allFiringRates','allMTPower','allFFTVals','allNumTrials','allTargetOnsetTimes','freqValsMT','freqValsFFT');
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Get Experimental Details %%%%%%%%%%%%%%%%%%%%%%%
    fileNameStringList = getAttentionExperimentDetails;
    fileNameStringListAll = cat(2,fileNameStringList{1},fileNameStringList{2});
    numSessions = length(fileNameStringListAll);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Get Timing Information %%%%%%%%%%%%%%%%%%%%%%%%%%
    timeRange = [-0.5 0];
    x = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{1}));
    timeVals = x.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    timePos = intersect(find(timeVals>=timeRange(1)),find(timeVals<timeRange(2)));
    numConditions = length(x.goodSpikeData);
    
    %%%%%%%%%%%%%%%%%%%%%% Set up MT and FFT parameters %%%%%%%%%%%%%%%%%%%%%%%
    fMax = 200;
    params.tapers = tapers;
    params.pad = -1;
    params.Fs = Fs;
    params.fpass = [0 fMax];
    params.trialave = 0;
    
    freqValsFFTall = 0:1/diff(timeRange):Fs-1/diff(timeRange);
    fPosFFT = (freqValsFFTall<=fMax);
    freqValsFFT = freqValsFFTall(fPosFFT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allFiringRates = cell(1,numSessions);
    allMTPower = cell(1,numSessions);
    allFFTVals = cell(1,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    
    for s=1:numSessions
        % Get Neural data
        data = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{s}));
        disp([num2str(s) ': ' fileNameStringListAll{s}]);
        
        clear firingRate mtPower fftVals
        firingRate = cell(2,numConditions);
        mtPower = cell(2,numConditions);
        fftVals = cell(2,numConditions);
        
        for c=1:numConditions
            for array=1:2
                numElectrodes = size(data.goodLFPData{array,c},1);
                
                for e=1:numElectrodes
                    firingRate{array,c}(e,:) = getSpikeCounts(data.goodSpikeData{array,c}(e,:),timeRange)./diff(timeRange);
                    
                    tmpLFP = squeeze(data.goodLFPData{array,c}(e,:,timePos))';
                    [mtPower{array,c}(e,:,:),freqValsMT] = mtspectrumc(tmpLFP,params);
                    fftx = fft(tmpLFP);
                    fftVals{array,c}(e,:,:) = fftx(fPosFFT,:);
                end
            end
        end
        allFiringRates{s} = firingRate;
        allMTPower{s} = mtPower;
        allFFTVals{s} = fftVals;
        allNumTrials{s} = data.nTrials;
        allTargetOnsetTimes{s} = data.targetOnsetTimes;
    end
    
    save(fileNameSave,'allFiringRates','allMTPower','allFFTVals','allNumTrials','allTargetOnsetTimes','freqValsMT','freqValsFFT');
end
end