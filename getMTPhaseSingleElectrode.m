% This program saves phases computed using MT

function [allFiringRates,allMTPhase,allNumTrials,allTargetOnsetTimes,freqValsMT] = getMTPhaseSingleElectrode(TWNum)

if ~exist('TWNum','var');                   TWNum=3;                    end

tapers = [TWNum 2*TWNum-1];

folderSavedData = fullfile(pwd,'savedData');
fileNameSaveMT = fullfile(folderSavedData,['singleElectrodeMTPhase' num2str(TWNum) '.mat']);

if exist(fileNameSaveMT,'file')

    disp(['Loading saved data in ' fileNameSaveMT]);
    load(fileNameSaveMT,'allFiringRates','allMTPhase','allNumTrials','allTargetOnsetTimes','freqValsMT');
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Set up MT  parameters %%%%%%%%%%%%%%%%%%%%%%%
    fMax = 200;
    params.tapers = tapers;
    params.pad = -1;
    params.Fs = Fs;
    params.fpass = [0 fMax];
    params.trialave = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allFiringRates = cell(1,numSessions);
    allMTPhase = cell(1,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    
    for s=1:numSessions
        % Get Neural data
        data = load(fullfile(folderSavedData,'neuralData',fileNameStringListAll{s}));
        disp([num2str(s) ': ' fileNameStringListAll{s}]);
        
        clear firingRate mtPhase
        firingRate = cell(2,numConditions);
        mtPhase = cell(2,numConditions);
        
        for c=1:numConditions
            for array=1:2
                numElectrodes = size(data.goodLFPData{array,c},1);
                
                for e=1:numElectrodes
                    firingRate{array,c}(e,:) = getSpikeCounts(data.goodSpikeData{array,c}(e,:),timeRange)./diff(timeRange);
                    
                    tmpLFP = squeeze(data.goodLFPData{array,c}(e,:,timePos))';
                    [~,freqValsMT,J] = mtspectrumc_returnJ(tmpLFP,params);
                    mtPhase{array,c}(e,:,:,:) = angle(J);
                end
            end
        end
        allFiringRates{s} = firingRate;
        allMTPhase{s} = mtPhase;
        allNumTrials{s} = data.nTrials;
        allTargetOnsetTimes{s} = data.targetOnsetTimes;
    end
    
    save(fileNameSaveMT,'allFiringRates','allMTPhase','allNumTrials','allTargetOnsetTimes','freqValsMT','-v7.3');
end
end