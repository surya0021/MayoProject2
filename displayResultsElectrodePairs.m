% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis
% measure - phase or power

function displayResultsElectrodePairs(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,measure)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('measure','var');                 measure='phase';            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the original conditions - AttLeftHit AttRightHit AttLeftMiss AttRightMiss
% colorNamesList(1,:) = [0 0 1]; % Attend-In Hit Blue
% colorNamesList(2,:) = [1 0 0]; % Attend-Out Hit Red
% colorNamesList(3,:) = [0 1 1]; % Attend-In Miss Cyan
% colorNamesList(4,:) = [1 0 1]; % Attend-Out Miss: Magenta

%%%%%%%%%%%%%%%%%% To plot differences and dPrimes %%%%%%%%%%%%%%%%%%%%%%%%
indicesToCompare{1} = [1 2]; % L vs R for Hits
indicesToCompare{2} = [3 4]; % L vs R for Misses
indicesToCompare{3} = [1 3]; % H vs M for Attend
indicesToCompare{4} = [2 4]; % H vs M for Attend R

colorsForComparison{1} = [0 1 0]; % Green
colorsForComparison{2} = [1 1 0]; % Yellow
colorsForComparison{3} = [0 0 1]; % Blue
colorsForComparison{4} = [1 0 0]; % Red

legendForComparison{1} = 'H(In-Out)';
legendForComparison{2} = 'M(In-Out)';
legendForComparison{3} = 'In(H-M)';
legendForComparison{4} = 'Out(H-M)';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlots = getPlotHandles(2,2,[0.05 0.05 0.9 0.9],0.05,0.1,0);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Check for intermediate data %%%%%%%%%%%%%%%%%%%
folderSavedData = fullfile(pwd,'savedData');
pairwiseDataToSave = fullfile(folderSavedData,['pairwiseData' conditionType num2str(targetOnsetMatchingChoice) 'N' num2str(numTrialCutoff) 'TW' num2str(TWNum) measure '.mat']);

if exist(pairwiseDataToSave,'file')
    disp(['Loading saved data in ' pairwiseDataToSave]);
    load(pairwiseDataToSave,'pairwiseMeasureData','freqValsMT');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(measure,'phase')
        [~,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getMTPhaseSingleElectrode(TWNum);
    elseif strcmp(measure,'power')
        [~,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getMTPowerSingleElectrode(TWNum);
    end
    numSessions = length(allTargetOnsetTimes0);
    numConditions = length(allTargetOnsetTimes0{1});
    
    %%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    targetTimeBinWidthMS = 250;
    goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);
    
    %%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
    allMTMeasure = cell(1,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    
    for i=1:numSessions
        tmpMTMeasure = cell(2,numConditions);
        tmpAllNumTrials = zeros(1,numConditions);
        tmpTargetOnsetTimes = cell(1,numConditions);
        
        for k=1:numConditions
            for j=1:2
                tmpMTMeasure{j,k} = allMTMeasure0{i}{j,k}(:,:,:,goodStimNums{i}{k});
            end
            tmpAllNumTrials(k) = length(goodStimNums{i}{k});
            tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
        end
        allMTMeasure{i} = tmpMTMeasure;
        allNumTrials{i} = tmpAllNumTrials;
        allTargetOnsetTimes{i} = tmpTargetOnsetTimes;
    end
    
    clear allMTMeasure0
    %%%%%%%%%%%%%%%%%%%%%%%%%% Get Condition Indices %%%%%%%%%%%%%%%%%%%%%%%%%%
    % The order of the 12 conditions is as follows: {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'};
    % Side 0 is Left (L), Side 1 is Right (R)
    fullOriginalConditionsList ={'HLV','HRV','HLI','HRI','MLV','MRV','MLI','MRI','HLN','HRN','MLN','MRN'};
    if strcmpi(conditionType(1),'V')
        conditionsToUse = [1 2 5 6];
    elseif strcmpi(conditionType(1),'N')
        conditionsToUse = 9:12;
    elseif strcmpi(conditionType(1),'I')
        conditionsToUse = [3 4 7 8];
    end
    originalConditionsList = fullOriginalConditionsList(conditionsToUse);
    numConditionsToUse = length(conditionsToUse);
    
    %%%%%%%%%%%%%%%%%%%% Select sessions with numTrials>=cutoff %%%%%%%%%%%%%%%
    allNumTrialsMatrix = cell2mat(allNumTrials');
    minTrialsConditions = min(allNumTrialsMatrix(:,conditionsToUse),[],2);
    
    badSessionList = find(minTrialsConditions<=numTrialCutoff);
    goodSessionList = setdiff(1:length(allNumTrials),badSessionList);
    numGoodSessions = length(goodSessionList);
    disp(['Discarded sessions: ' num2str(badSessionList')]);

%     %%%%%%%%%%%%%%%%%%%%%%% Display TargetOnset Histogram %%%%%%%%%%%%%%%%%%%%%
%     targetOnsetTimesForHistogram = cell(1,numConditionsToUse);
%     for i=1:numGoodSessions
%         for j=1:numConditionsToUse
%             targetOnsetTimesForHistogram{j} = cat(2,targetOnsetTimesForHistogram{j},allTargetOnsetTimes{goodSessionList(i)}{conditionsToUse(j)});
%         end
%     end
%     
%     targetOnsetEdges = 500:targetTimeBinWidthMS:5500;
%     c = 500+targetTimeBinWidthMS/2:targetTimeBinWidthMS:5500;
%     legendStr = cell(1,numConditionsToUse);
%     for i=1:numConditionsToUse
%         h = histcounts(targetOnsetTimesForHistogram{i},targetOnsetEdges);
%         plot(hPlots(1,1),c,h,'color',colorNamesList(i,:)); hold(hPlots(1,1),'on');
%         legendStr{i} = [originalConditionsList{i} '(' num2str(length(targetOnsetTimesForHistogram{i})) ')'];
%     end
%     legend(hPlots(1,1),legendStr);
%     xlabel(hPlots(1,1),'TargetOnset (ms)'); ylabel(hPlots(1,1),'Num Stim');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    goodMTMeasure = allMTMeasure(goodSessionList);
    pairwiseMeasureData = cell(2,numConditionsToUse,numGoodSessions);
    for side=1:2
        for c=1:numConditionsToUse
            for i=1:numGoodSessions
                disp(['side' num2str(side) '_' originalConditionsList{c} '_session' num2str(i)]);
                tmpData = goodMTMeasure{i}{side,conditionsToUse(c)};
                numElecs = size(tmpData,1);
                
                tmpMeasureList = [];
                for e1=1:numElecs
                    for e2=e1+1:numElecs
                        d1 = squeeze(tmpData(e1,:,:,:));
                        d2 = squeeze(tmpData(e2,:,:,:));
                        tmpMeasureList = cat(3,tmpMeasureList,getMeasure(d1,d2,measure));
                    end
                end
                pairwiseMeasureData{side,c,i} = tmpMeasureList;
            end
        end
    end
    save(pairwiseDataToSave,'pairwiseMeasureData','freqValsMT');
end

%%%%%%%%%%%%%%%%%%%% Get Means and DPrimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pairwise data is for conditions H0, H1, M0 and M1. This needs to be
% converted to Hit(InVsOut), Miss(InVsOut), In(HitVsMiss) and Out(HitVsMiss)
allMeans0 = cell(2,4);
allDPrimes0 = cell(2,4);
numGoodSessions = size(pairwiseMeasureData,3);
for side=1:2
    for c=1:4 % Comparison
        disp(['side' num2str(side) '_' legendForComparison{c}]);      
        tmpIndexToCompare = indicesToCompare{c};
    
        meanDataAllSessions=[];
        dPrimeDataAllSessions=[];
        for i=1:numGoodSessions
            tmpData1 = pairwiseMeasureData{side,tmpIndexToCompare(1),i};
            tmpData2 = pairwiseMeasureData{side,tmpIndexToCompare(2),i};
            
            [numFreqs,~,numPairs] = size(tmpData1);
            
            clear tmpMeanData tmpDPrime
            tmpMeanData = zeros(numFreqs,numPairs);
            tmpDPrimeData = zeros(numFreqs,numPairs);
            for f=1:numFreqs
                for p=1:numPairs
                    x1 = squeeze(tmpData1(f,:,p));
                    x2 = squeeze(tmpData2(f,:,p));
                    tmpMeanData(f,p) = mean(x1)-mean(x2);
                    tmpDPrimeData(f,p) = getDPrime(x1,x2);
                end
            end
            meanDataAllSessions = cat(2,meanDataAllSessions,tmpMeanData);
            dPrimeDataAllSessions = cat(2,dPrimeDataAllSessions,tmpDPrimeData);
        end
        allMeans0{side,c} = meanDataAllSessions;
        allDPrimes0{side,c} = dPrimeDataAllSessions;
    end
end

allMeans = combineComparisonDataAcrossBothArrays(allMeans0);
allDPrimes = combineComparisonDataAcrossBothArrays(allDPrimes0);

% Plot DPrimes
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    plotData(hPlots(plotPos,1),freqValsMT,allMeans{i}',colorsForComparison{i},1);
    plotData(hPlots(plotPos,2),freqValsMT,allDPrimes{i}',colorsForComparison{i},1);
end
for i=1:2
    plot(hPlots(i,2),freqValsMT,zeros(1,length(freqValsMT)),'k--');
    if strcmp(measure,'phase')
        ylabel(hPlots(i,1),'\Delta PPC');
    else
        ylabel(hPlots(i,1),'\Delta Corr');
    end
    plot(hPlots(i,2),freqValsMT,zeros(1,length(freqValsMT)),'k--');
    ylabel(hPlots(i,2),'dPrime');
end
end

function outMeasure = getMeasure(d1,d2,measure)

[numFreqs,~,numTrials] = size(d1);

outMeasure = zeros(numFreqs,numTrials);
for i=1:numFreqs
    for t=1:numTrials
        tmpData1 = squeeze(d1(i,:,t));
        tmpData2 = squeeze(d2(i,:,t));
        if strcmp(measure,'phase') % Get PPC
            outMeasure(i,t) = getPPC(tmpData1-tmpData2);
        elseif strcmp(measure,'power') % Get Correlation
        end
    end
end
end
function combinedData = combineComparisonDataAcrossBothArrays(data)

combinedData = cell(1,4);

% Input data has size of 2x4. Rows correspond to sides L and R.
% The 4 comparisons are 1-2, 3-4, 1-3 and 3-4, corresponding to H0-H1,
% M0-M1, H0-M0 and H1-M1. Here 0 and 1 correspond to sides L and R.

% AttIn-AttOut Hit & Miss
% The comparison 1-2 (H0-H1) corresponds to AttIn-out Hit for side 0 and
% AttOut - AttIn Hit for side 1. Therefore We can get AttIn-Out Hit by
% simply flipping the values for side 1 and concatenating

combinedData{1} = [data{1,1} -data{2,1}];
combinedData{2} = [data{1,2} -data{2,2}];

% Hits-Miss for AttIn and AttOut: 
% For comparison 3 (1-3 or H0-M0) - these correspond to AttInHvsM for array
% 0 but AttOutHVsM for array 1. This flips for comparison 4.

combinedData{3} = [data{1,3} data{2,4}];
combinedData{4} = [data{1,4} data{2,3}];

end
function plotData(hPlot,xs,data,colorName,showSEMFlag)

if ~exist('showSEMFlag','var');     showSEMFlag = 1;                    end

tmp=rgb2hsv(colorName); 
tmp2 = [tmp(1) tmp(2)/3 tmp(3)];
colorName2 = hsv2rgb(tmp2); % Same color with more saturation

mData = squeeze(mean(data,1));

if showSEMFlag 
    sData = std(data,[],1)/sqrt(size(data,1));    
    xsLong = [xs fliplr(xs)];
    ysLong = [mData+sData fliplr(mData-sData)];
    patch(xsLong,ysLong,colorName2,'EdgeColor','none','parent',hPlot);
end
hold(hPlot,'on');
plot(hPlot,xs,mData,'color',colorName,'linewidth',2); 
end