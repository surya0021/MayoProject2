% This function computes single trial measures of shuffle corrected and raw (uncorrected) power correlation 
% and PPC using hilbert transform

% Input:
% freqHalfBandwidth - Half the frequency bandwidth of LFP butterworth filter 

function [ppcHT,powerCorrHT,timeValsHT,freqValsHT,freqHalfBandwidth,allNumTrials,allTargetOnsetTimes] = getHilbertPairwiseMeasures(freqHalfBandwidth)
if ~exist('freqHalfBandwidth','var');          freqHalfBandwidth = 5;   end
folderSavedData = fullfile(pwd,'savedData');

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Experimental Details %%%%%%%%%%%%%%%%%%%%%%%
fileNameStringList = getAttentionExperimentDetails;
fileNameStringListAll = cat(2,fileNameStringList{1},fileNameStringList{2});
numSessions = length(fileNameStringListAll);

%%%%%%%%%%%%%%%%%%%%%%%%% Get Timing Information %%%%%%%%%%%%%%%%%%%%%%%%%%
timeRange = [-0.5 0];
x = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{1}));
timeVals = x.timeVals;
timePos = intersect(find(timeVals>=timeRange(1)),find(timeVals<timeRange(2)));
numConditions = length(x.goodSpikeData);
fullOriginalConditionsList ={'HLV','HRV','HLI','HRI','MLV','MRV','MLI','MRI','HLN','HRN','MLN','MRN'};
% HT parameters
timeValsHT = timeVals(timePos);
Fs = round(1/(timeValsHT(2)-timeValsHT(1)));
fPass = [0 200];
% freqValsHT = [freqHalfBandwidth+1 2*freqHalfBandwidth:freqHalfBandwidth:fPass(2)];
%freqValsHT = [freqHalfBandwidth+1 freqHalfBandwidth+2:2:fPass(2)]; % Overlapping for even freqHalfBandwidth
freqValsHT = freqHalfBandwidth+1:2:fPass(2); % Overlapping for odd freqHalfBandwidth
filterOrder = 4;

% fileNameSave = [];
fileNameSave = fullfile(folderSavedData,['ElectrodePairwiseMeasuresHTW' num2str(freqHalfBandwidth) '.mat']);
if exist(fileNameSave,'file')
    load(fileNameSave) %#ok<LOAD> 
else
    ppcHT = cell(2,numConditions,numSessions);
    powerCorrHT = cell(2,numConditions,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    for s=1:numSessions
        % Get Neural data
        data = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{s}));
        disp([num2str(s) ': ' fileNameStringListAll{s}]);

        clear htPower htPhase
        htVals = cell(2,numConditions);
        for c=1:numConditions
            disp(['Condition: ' fullOriginalConditionsList{c}])
            for side=1:2
                numEles = size(data.goodLFPData{side,c},1);
                for e=1:numEles
                    tmpLFP = squeeze(data.goodLFPData{side,c}(e,:,timePos))';
                    htVals{side,c}(:,:,:,e) = getHilbert(tmpLFP,Fs,freqValsHT,freqHalfBandwidth,filterOrder);
                end
                ppcTMP = [];
                powerCorrTMP = [];
                for e=1:numEles
                    for f=e+1:numEles
                        phase1 = angle(htVals{side,c}(:,:,:,e));
                        phase2 = angle(htVals{side,c}(:,:,:,f));
                        power1 = abs(htVals{side,c}(:,:,:,e)).^2;
                        power2 = abs(htVals{side,c}(:,:,:,f)).^2;
                        ppcTMP = cat(3,ppcTMP,getMeasures(phase1,phase2,'phase'));
                        powerCorrTMP = cat(3,powerCorrTMP,getMeasures(power1,power2,'power'));
                    end
                end
                ppcHT{side,c,s} = ppcTMP;
                powerCorrHT{side,c,s} = powerCorrTMP;
            end
        end
        allNumTrials{s} = data.nTrials;
        allTargetOnsetTimes{s} = data.targetOnsetTimes;
    end
    save(fileNameSave,'ppcHT','powerCorrHT','timeValsHT','freqValsHT','freqHalfBandwidth','allNumTrials','allTargetOnsetTimes','-v7.3');
end

end

function htVals = getHilbert(data,Fs,freqVals,freqHalfBandwidth,filterOrder)

for f=1:length(freqVals)
    clear bandPassSignal
    freqRange = [freqVals(f)-freqHalfBandwidth freqVals(f)+freqHalfBandwidth];
    [b,a] = butter(filterOrder,freqRange/(Fs/2),'bandPass');
    bandPassSignal = filtfilt(b,a,data);
    htVals(:,f,:) = hilbert(bandPassSignal); %#ok<AGROW> 
end
% htPower = (abs(htVals)).^2;
% htPhase = angle(htVals);
end

function outMeasure = getMeasures(data1,data2,measure)
nTrials = size(data1,3);
nFreqs = size(data1,2);

if nTrials~=size(data2,3)
    error('Trials does not match')
end

if strcmp(measure,'phase')
    outMeasure = zeros(nTrials,nFreqs);
    for f=1:nFreqs
        for t=1:nTrials
            outMeasure(t,f) = getPPC(data1(:,f,t)-data2(:,f,t));
        end
    end
elseif strcmp(measure,'power')
    outMeasure = zeros(nTrials,nFreqs,1,3); % Last dimension (3) is for corrected,raw and shuffled correlation values respectively
    for f=1:nFreqs
        outMeasureTMP = corr(squeeze(data1(:,f,:)),squeeze(data2(:,f,:))); %
        outMeasure(:,f,1,2) = diag(outMeasureTMP); % Raw correlations
        ind = 1:size(outMeasureTMP,1)+1:size(outMeasureTMP,1)^2; % indices for diagonal elements
        if isequaln(outMeasureTMP(ind)',outMeasure(:,f,1,2))
            outMeasureTMP(ind)=nan;
        else
            error('Mismatch between diagonal elements and its index')
        end
        outMeasure(:,f,1,3) = mean(outMeasureTMP,2,'omitnan'); % shuffle correlation across non-simultaneous trials 
        outMeasure(:,f,1,1) = outMeasure(:,f,1,2)-outMeasure(:,f,1,3); % Corrected correlation
    end
end
end
