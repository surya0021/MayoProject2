% This program follows the same logic as getMTValsSingleElectrode, but with important differences. 
% When using MT, power and phase from each taper is saved separately,
% yielding multiple values per trial, which are later used for single trial
% analysis. Here, each trial is first divided into numDivisions segments.
% Then, firing rate and FFT power and phase is computed for each segment.
% Note that the frequency resolution will decrease by numDivisions.

function [allFiringRates,allFFTVals,allNumTrials,allTargetOnsetTimes,freqValsFFT] = getFFTValsSingleElectrode(numDivisions,measure)

if ~exist('numDivisions','var');           numDivisions=10;             end

folderSavedData = fullfile(pwd,'savedData');

if strcmp(measure,'phase')
    fileNameSaveFFT = fullfile(folderSavedData,['singleElectrodeFFTPhase' num2str(numDivisions) '.mat']);
else
    fileNameSaveFFT = fullfile(folderSavedData,['singleElectrodeFFTPower' num2str(numDivisions) '.mat']);
end

if exist(fileNameSaveFFT,'file')

    disp(['Loading saved data in ' fileNameSaveFFT]);
    
    if strcmp(measure,'phase')
        load(fileNameSaveFFT,'allFiringRates','allFFTPhase','allNumTrials','allTargetOnsetTimes','freqValsFFT');
        allFFTVals = allFFTPhase;
    elseif strcmp(measure,'power')
        load(fileNameSaveFFT,'allFiringRates','allFFTPower','allNumTrials','allTargetOnsetTimes','freqValsFFT');
        allFFTVals = allFFTPower;
    end
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Get Experimental Details %%%%%%%%%%%%%%%%%%%%%%%
    fileNameStringList = getAttentionExperimentDetails;
    fileNameStringListAll = cat(2,fileNameStringList{1},fileNameStringList{2});
    numSessions = length(fileNameStringListAll);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Get Timing Information %%%%%%%%%%%%%%%%%%%%%%%%%%
    timeRange = [-0.5 0]; % full time range
    segmentDur = diff(timeRange)/numDivisions;
    segmentRangeList = cell(1,numDivisions);
    for i=1:numDivisions
        segmentRangeList{i} = (timeRange(1) + (i-1)*segmentDur) + [0 segmentDur];
    end
    
    x = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{1}));
    timeVals = x.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    segmentLength = round(segmentDur*Fs);
    
    numConditions = length(x.goodSpikeData);
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Set up FFT  parameters %%%%%%%%%%%%%%%%%%%%%%
    fMax = 200;
    freqValsFFTFull = 0:(1/segmentDur):Fs-(1/segmentDur);
    goodFPos = 1:find(freqValsFFTFull>=fMax,1);
    freqValsFFT = freqValsFFTFull(goodFPos);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allFiringRates = cell(1,numSessions);
    allFFTPhase = cell(1,numSessions);
    allFFTPower = cell(1,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    
    for s=1:numSessions
        % Get Neural data
        data = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{s}));
        disp([num2str(s) ': ' fileNameStringListAll{s}]);
        
        clear firingRate fftPhase fftPower
        firingRate = cell(2,numConditions);
        fftPhase = cell(2,numConditions);
        fftPower = cell(2,numConditions);
        
        for c=1:numConditions
            for array=1:2
                numElectrodes = size(data.goodLFPData{array,c},1);            
                for e=1:numElectrodes
                    for n=1:numDivisions
                        thisSegmentRange = segmentRangeList{n};
                        timePos = find(timeVals<=thisSegmentRange(1),1,'last') + (0:segmentLength-1);
                        
                        firingRate{array,c}(e,1,n,:) = getSpikeCounts(data.goodSpikeData{array,c}(e,:),thisSegmentRange)./diff(thisSegmentRange);
                    
                        tmpLFP = squeeze(data.goodLFPData{array,c}(e,:,timePos))';
                        tmpFFT = fft(tmpLFP); tmpFFT = tmpFFT(goodFPos,:);
                        fftPhase{array,c}(e,:,n,:) = angle(tmpFFT);
                        fftPower{array,c}(e,:,n,:) = abs(tmpFFT).^2;
                    end
                end
            end
        end
        allFiringRates{s} = firingRate;
        allFFTPhase{s} = fftPhase;
        allFFTPower{s} = fftPower;
        allNumTrials{s} = data.nTrials;
        allTargetOnsetTimes{s} = data.targetOnsetTimes;
    end
    
    if strcmp(measure,'phase')
        save(fileNameSaveFFT,'allFiringRates','allFFTPhase','allNumTrials','allTargetOnsetTimes','freqValsFFT','-v7.3');
        allFFTVals=allFFTPhase;
    elseif strcmp(measure,'power')
        save(fileNameSaveFFT,'allFiringRates','allFFTPower','allNumTrials','allTargetOnsetTimes','freqValsFFT','-v7.3');
        allFFTVals=allFFTPower;
    end
end
end