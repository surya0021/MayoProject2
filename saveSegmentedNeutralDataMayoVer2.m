% This program is second version of segmentData program used for
% decoding project where the the time was 500 ms before target onset. Here
% we segment the data to 1000 ms before target onset
% SP 27/10/2021


%This function saves the neutral trials which are divided further depending
%on target location
%Here the format is {trialOutcome}{targetLocation}{cueType} used instead of
%previous format {trialOutcome}{cueLocation}{cueType} c.f.
%savesegmentedDataMayo.m
%{trialOutcome}= Hit(H) or Miss(M);  {targetLocation}= Left(0) or Right(1);
%{cueType}= Neutral(N)

function saveSegmentedNeutralDataMayoVer2(fileNameString,folderSourceString)

if ~exist('folderSourceString','var');   folderSourceString='E:\Mayo';       end
neutralTrialFlag=1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
saveStringConditionList = [{'H0N'} {'H1N'} {'M0N'} {'M1N'}]; % Data must be stored in this order
Fs = 2000;
timePeriodS = [-1 0.1]; % Save this interval around the target onset
saveStringTimePeriodList = {'_TargetOnset'};

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
goodIndexList = getGoodIndices(CDS,DAT,[],neutralTrialFlag); % Gets Indices for the 12 categories
trialStartTimeS = cell2mat(cellfun(@(x) x{1}, CDS(:,1),'UniformOutput',false ));
stimulusOnTimeS = cellfun(@(x) x{5}, CDS(:,1),'UniformOutput',false); % still under cell format. dimension: trials x 1 (in the first and only columnm there are cells containing all the trial start times for the corresponding trial
%saccadeTimeS = cell2mat(cellfun(@(x) x{7}, CDS(:,1),'UniformOutput',false));

[~,~,targetOnTimeMS,~,~] = getInfoDATFile(DAT);

numElectrodes = size(lfpData,2);
numTimePos = round(Fs*diff(timePeriodS));


for i=1:length(saveStringConditionList)
    disp(['Working on condition ' saveStringConditionList{i}]);
    
    
        clear segmentedLFPData segmentedSpikeData
                
        fileNameSaveLFP = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '_LFP.mat']);
        fileNameSaveSpikes = fullfile(folderNameOut,[fileNameString saveStringConditionList{i} saveStringTimePeriodList num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '_Spikes.mat']);
        fileNameSaveGoodIndList = fullfile(folderNameOut,['goodIndexList_SplitNeutral_TargetOnset' num2str(timePeriodS(1)) '_' num2str(timePeriodS(2)) '.mat']);
        
        numTrials = length(goodIndexList{8+i});  % GoodIndexList for Neutral Trials starts from the 9th cell 9:H0N , 10:H1N , 11:M0N, 12:M1N
%         segmentedLFPData = zeros(numElectrodes,numTrials,numTimePos(j));
%         segmentedSpikeData = cell(numElectrodes,numTrials);
        kPrime = 0;
        for k=1:numTrials
            t=goodIndexList{8+i}(k); % Trial Number of interest

            startPosLFP = round(Fs*(targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1)-trialStartTimeS(t) + timePeriodS(1)));
            startPosSpikes = targetOnTimeMS(t)/1000 + stimulusOnTimeS{t}(1);
            
            if targetOnTimeMS(t)>=1250     % Choosing the trials in which target occurs 1250 ms after stimulus. Additional 250 ms from 1000 ms is set to avoid the stimulus onset transient response.
               
                kPrime = kPrime+1;
                goodIndexListNew{i}(kPrime) = t; %#ok<AGROW>
                for n=1:numElectrodes
                    segmentedLFPData(n,kPrime,:) = lfpData{t,n}(startPosLFP+(1:numTimePos)); %#ok<AGROW>
                    spikeDataThisTrial = spikeData{n,1,t};
                    if ~isempty(spikeDataThisTrial)
                        spikeTimesThisTrial = spikeDataThisTrial - startPosSpikes; % times relative to stimulus or target onset
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