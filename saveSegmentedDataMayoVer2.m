% This program is second version of the segmentData program used for
% decoding project where the time was 500 ms before target onset. Here
% data is segmented 1000 ms before target onset
% SP 23/10/2021


% This program takes the fileNameString and saves the segmented data in
% appropriate folders

% fileNameString: string indicating the fileName.
% We assume that the following files are available at
% folderSourceString\Data\extractedData
% 1. fileNameString_Codes
% 2. fileNameString_DAT
% 3. fileNameString_extractedTrialsLFP

% Followng Patrick's nomenclature, we save data in the format
% {trialOutcome}{CueType}{CueLocation}.
% {trialOutcome}: Hit (H) or Miss (M)
% {CueType}: Valid (V), Invalid (I) or Neutral (N)
% {CueLocation}: In case the CueType is Valid or Invalid, we also have the
% cue location (0: left or 1: right).
% So there are the following 10 file types: H0V,H1V,H0I,H1I,M0V,M1V,M0I,M1I,HN,MN

% Further, we save data around target.

function saveSegmentedDataMayoVer2(fileNameString,folderSourceString,instructionTrialFlag)

if ~exist('folderSourceString','var');   folderSourceString = 'E:\Mayo';   end     %'C:\Supratim\Projects\MayoProject\';    
if ~exist('instructionTrialFlag','var');    instructionTrialFlag=0;     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
saveStringConditionList = [{'H0V'} {'H1V'} {'H0I'} {'H1I'} {'M0V'} {'M1V'} {'M0I'} {'M1I'} {'HN'} {'MN'}]; % Data must be stored in this order
Fs = 2000;
timePeriodS = [-1 0.1]; % Save this interval around the target onset
saveStringTimePeriodList =  {'_TargetOnset'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderNameIn = fullfile(folderSourceString,'Data','extractedData');
folderNameOut = fullfile(folderSourceString,'Data','segmentedData2',fileNameString);
makeDirectory(folderNameOut);

fileNameCDS = [fileNameString '_Codes'];
CDS = load(fullfile(folderNameIn,fileNameCDS));
CDS = CDS.(fileNameCDS);
disp([fileNameCDS ' loaded....']);

fileNameDAT = strcat(fileNameString, '_DAT');
DAT = load(fullfile(folderNameIn,fileNameDAT));
DAT = DAT.(fileNameDAT);
disp([fileNameDAT ' loaded....']);

fileNameLFP = strcat(fileNameString, '_extractedTrialsLFP');
lfpData = load(fullfile(folderNameIn,fileNameLFP));
if isfield(lfpData,fileNameLFP)
    lfpData = lfpData.(fileNameLFP); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
else
    lfpData = lfpData.extractedTrialsLFP;
end
disp([fileNameLFP ' file loaded.....']);

fileNameSpikes = strcat(fileNameString, '_Spikes');
spikeData = load(fullfile(folderNameIn,fileNameSpikes));
spikeData = spikeData.(fileNameSpikes); % dimension: trials x channels (each cell is an array of the voltage values for the specific trial and channel)
disp([fileNameSpikes ' file loaded.....']);

%%%%%%%%%%%%%%%%%%%%%%%% Find appropriate indices %%%%%%%%%%%%%%%%%%%%%%%%%
goodIndexList = getGoodIndices(CDS,DAT,instructionTrialFlag); % Get Indices for the 10 categories
trialStartTimeS = cell2mat(cellfun(@(x) x{1}, CDS(:,1),'UniformOutput',false ));
stimulusOnTimeS = cellfun(@(x) x{5}, CDS(:,1),'UniformOutput',false); % still under cell format. dimension: trials x 1 (in the first and only columnm there are cells containing all the trial start times for the corresponding trial
%saccadeTimeS = cell2mat(cellfun(@(x) x{7}, CDS(:,1),'UniformOutput',false));

[~,~,targetOnTimeMS,~,~] = getInfoDATFile(DAT);

numElectrodes = size(lfpData,2);
numTimePos = round(Fs*diff(timePeriodS));


for i=1:length(goodIndexList) % For each of the 10 conditions
    disp(['Working on condition ' saveStringConditionList{i}]);
    
    clear segmentedLFPData segmentedSpikeData
    
    if ~instructionTrialFlag
        fileNameSaveLFP = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '_LFP.mat']);
        fileNameSaveSpikes = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '_Spikes.mat']);
        fileNameSaveGoodIndList = fullfile(folderNameOut,['goodIndexList_TargetOnset' num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '.mat']);
    else
        fileNameSaveLFP = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList '_LFP_Instruct.mat']);
        fileNameSaveSpikes = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList '_Spikes_Instruct.mat']);
    end
    numTrials = length(goodIndexList{i});
    %         segmentedLFPData = zeros(numElectrodes,numTrials,numTimePos(j));
    %         segmentedSpikeData = cell(numElectrodes,numTrials);
    
    kPrime = 0;     % Initializing an index for trial number
    for k=1:numTrials
        t=goodIndexList{i}(k); % Trial Number of interest
        if targetOnTimeMS(t)>=1250  % Choosing the trials in which target occurs 1250 ms after stimulus. Additional 250 ms from 1000 ms is set to avoid the stimulus onset transient response.
            kPrime = kPrime+1;
            goodIndexListNew{i}(kPrime) = t; %#ok<AGROW>
            startPosLFP = round(Fs*(targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1)-trialStartTimeS(t) + timePeriodS(1)));
            startPosSpikes = targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1);
            
            
            for n=1:numElectrodes
                segmentedLFPData(n,kPrime,:) = lfpData{t,n}(startPosLFP+(1:numTimePos)); %#ok<AGROW>
                spikeDataThisTrial = spikeData{n,1,t};
                if ~isempty(spikeDataThisTrial)
                    spikeTimesThisTrial = spikeDataThisTrial - startPosSpikes; % times relative target onset
                    usefulSpikePos = intersect(find(spikeTimesThisTrial>=timePeriodS(1)),find(spikeTimesThisTrial<timePeriodS(2)));
                    segmentedSpikeData{n,kPrime} = spikeTimesThisTrial(usefulSpikePos); %#ok<AGROW>
                end
            end
        end
    end
    timeVals = timePeriodS(1):1/Fs:timePeriodS(2)-1/Fs;
    save(fileNameSaveLFP,'timeVals','segmentedLFPData');
    save(fileNameSaveSpikes,'segmentedSpikeData');
end
save(fileNameSaveGoodIndList,'goodIndexListNew');
end