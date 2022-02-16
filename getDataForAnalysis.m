%This program saves LFP and spike data of a single orientation change chosen either among selected (the ones used in decoding project)
%or all orientation change by setting the "oriFlag" of all session and 12 attention condition in a single
%file which is used in all further analysis

%Prerequisite data : 
% 1)extractedData  2)segmentedData  3)electrodeArrayList

%Parameters to set 
% 1)OriFlag     2)targetOnTimeCriteria     3)targetOnTImeCutoff
folderSourceString = 'G:';
folderNameExtractedData = fullfile(folderSourceString,'Mayo','Data','extractedData');

oriFlag = 1; % 0: all orientation change are considered    1: only selected orientation change (which were used in decoding project) are considered.

targetOnTimeCriteria = 0; % 0 = No cutoff on target onset time  1: Cutoff on target onset time
targetOnTimeCutoff = 0; % cutoff in MS chosen for target onset time. Set this to 0 if targetOnTimeCriteria is set to 0

signalDuration = [-0.5 0.1]; % the time points between which the data is segmented.

fileNameStringList = getAttentionExperimentDetails;
fileNameStringListAll = cat(2,fileNameStringList{1},fileNameStringList{2});

newOrderForOriChangeIndex = [1:4,1:4,5:6,5:6];   % rearranging and repeating the Ori change index to match attCueList
newIndex = [1 2 9 10 5 6 11 12 3 4 7 8];    % Rearranging the attention conditions
attCueList = {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'};
newAttCueList = {'H0V','H1V','H0N','H1N','M0V','M1V','M0N','M1N','H0I','H1I','M0I','M1I'};

electrodeFolderString = fullfile(folderSourceString,'Mayo','Data','savedDataSummary');
electrodeFileString = fullfile(electrodeFolderString,'electrodeArrayListStimulated.mat');
eListAll = load(electrodeFileString);

saveDataFolder = fullfile(folderSourceString,'Projects','MayoProject2','Data','savedDataForAnalysis');
makeDirectory(saveDataFolder)
if oriFlag==0
    oriType = 'All';
elseif oriFlag==1
    oriType = 'Selected';
end
fileNameSaveString = ['data_' oriType 'Ori' '_Perform50' '_targetOnTimeCriterion' num2str(targetOnTimeCriteria) '_targetOnTimeCutoff' num2str(targetOnTimeCutoff) '_targetOnset' num2str(signalDuration(1)) '_' num2str(signalDuration(2)) '.mat']; 
fileNameSave = fullfile(saveDataFolder,fileNameSaveString);

goodLFPData = cell(1,2);
goodSpikeData = cell(1,2);
nTrials = zeros(length(fileNameStringListAll),length(attCueList));
for i=1:length(fileNameStringListAll)
    
    disp(['Session ' num2str(i) ' of ' num2str(length(fileNameStringListAll))])
    
    fileNameDAT = strcat(fileNameStringListAll{i}, '_DAT');
    DAT = load(fullfile(folderNameExtractedData,fileNameDAT));
    DAT = DAT.(fileNameDAT);
    [~,~,targetOnTimeMS,~,~] = getInfoDATFile(DAT);
    
    dataFolderString = fullfile(folderSourceString,'Mayo','Data','segmentedData',fileNameStringListAll{i});
    
    [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg] = getBehavior(fileNameStringListAll{i});
     
    if oriFlag == 0
        [~,oriChangeIndex] = min(abs(perCorrect-0.5),[],2);     % choosing orientation change near 50% performance
    elseif oriFlag == 1
        [~,oriChangeIndexTMP] = min(abs(perCorrect(:,[2 3])-0.5),[],2);
        oriChangeIndex = oriChangeIndexTMP+1;   % One is added to get the correct index for orientation change
    end
    
    allConditionOriChangeIndex = oriChangeIndex(newOrderForOriChangeIndex);
    eList = eListAll.electrodeArrayList{i};
    
    clear goodLFPDataTMP goodSpikeDataTMP
    goodLFPDataTMP = cell(1,2);
    goodSpikeDataTMP = cell(1,2);
    for j=1:length(attCueList)
        fileNameStringLFP = fullfile(dataFolderString,[fileNameStringListAll{i} attCueList{j} '_TargetOnset_' 'LFP.mat']);
        fileNameStringSpikes = fullfile(dataFolderString,[fileNameStringListAll{i} attCueList{j} '_TargetOnset_' 'Spikes.mat']);
        
        lfpData = load(fileNameStringLFP);
        spikeData = load(fileNameStringSpikes);
        
        oriChangeThisCondition = orientationChangeDeg(goodIndexList{j});
        targetOnTimeThisCondition = targetOnTimeMS(goodIndexList{j}); 
        goodPosTMP = find(oriChangeThisCondition==uniqueOrientationChangeDeg(allConditionOriChangeIndex(j)));
        
        if targetOnTimeCriteria == 1
            goodTargetPos = find(targetOnTimeThisCondition>=targetOnTimeCutoff);
            goodPos = intersect(goodPosTMP,goodTargetPos);
        elseif targetOnTimeCriteria == 0
            goodPos = goodPosTMP;
        end
        
        for array=1:2   % 1: Right Array    2: Left Array
            goodLFPDataTMP{array}{j} = lfpData.segmentedLFPData(eList{array},goodPos,:);
            goodSpikeDataTMP{array}{j} = spikeData.segmentedSpikeData(eList{array},goodPos);
        end
    end  
    %Rearranging the data array 
    for array=1:2
        goodLFPData{array}(i,:) = goodLFPDataTMP{array}(newIndex);
        goodSpikeData{array}(i,:) = goodSpikeDataTMP{array}(newIndex);
    end
    
    nTrials(i,:) = cellfun(@(x) size(x,2),goodLFPData{1}(i,:),'un',1);
    timeVals = lfpData.timeVals;
end

save(fileNameSave,'goodLFPData','goodSpikeData','timeVals','nTrials','-v7.3')
