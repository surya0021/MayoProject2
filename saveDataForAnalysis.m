% This program saves LFP and spike data for a single orientation change.

% folderSourceString - parent folder. % Data is assumed to be stored in {folderSourceString}/Data
% goodPos - the orientation change index to be saved. This should a 1x6 array, with indices for the 6 attention conditions (0V, 1V, 0I, 1I, 0N and 1N)
% electrodeList - list of useful electrodes for each side

function saveDataForAnalysis(fileNameString,folderSourceString,goodOriPos,electrodeList)

dataFolderString = fullfile(folderSourceString,'Data','segmentedData',fileNameString);
folderSave = fullfile(pwd,'savedData','neuralData');
makeDirectory(folderSave);
fileNameSave = fullfile(folderSave,fileNameString); % Output is saved in a local folder

% We use 6 conditions - valid cued attend side 0 or side 1 (0V and 1V),
% invalid cues for the same (0I and 1I), and neutral cue with target on
% side 0 or 1 (0N and 1N). For each, we have both hit and miss conditions.

[~,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg,~,targetOnTimeMS] = getBehavior(fileNameString);
attCueList = {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'}; % this is the order in which getBehavior returns goodIndexList 
allConditionsOriChangeIndex = [goodOriPos(1:4) goodOriPos(1:4) goodOriPos(5:6) goodOriPos(5:6)]; % Selected orientation in each of the 12 conditions - hits and misses have the same ori change

numConditions = length(attCueList);
goodLFPData = cell(2,numConditions);
goodSpikeData = cell(2,numConditions);
nTrials = zeros(1,numConditions);
targetOnsetTimes = cell(1,numConditions);

for i=1:numConditions
    
    % Get data from the main data folder
    fileNameStringLFP = fullfile(dataFolderString,[fileNameString attCueList{i} '_TargetOnset_LFP.mat']);
    fileNameStringSpikes = fullfile(dataFolderString,[fileNameString attCueList{i} '_TargetOnset_Spikes.mat']);
    lfpData = load(fileNameStringLFP);
    spikeData = load(fileNameStringSpikes);
    
    % Find orientation change for all stimulus repeats for this condition
    oriChangeThisCondition = orientationChangeDeg(goodIndexList{i});
    targetOnTimeThisCondition = targetOnTimeMS(goodIndexList{i});
    goodOriPos = find(oriChangeThisCondition==uniqueOrientationChangeDeg(allConditionsOriChangeIndex(i)));
    nTrials(i) = length(goodOriPos);
    targetOnsetTimes{i} = targetOnTimeThisCondition(goodOriPos);
    
    for array=1:2   % 1: Right Array    2: Left Array
        goodLFPData{array}{i} = lfpData.segmentedLFPData(electrodeList{array},goodOriPos,:);
        goodSpikeData{array}{i} = spikeData.segmentedSpikeData(electrodeList{array},goodOriPos);
    end
end

timeVals = lfpData.timeVals;
save(fileNameSave,'goodLFPData','goodSpikeData','timeVals','targetOnsetTimes','nTrials');
end
