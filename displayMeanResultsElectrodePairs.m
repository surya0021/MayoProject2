% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis
% measure - phase or power

% Modified from displayResultsElectrodePairs. Here, we compute PPC or
% correlation across trials also (yielding measures similar to previous paper)
function displayMeanResultsElectrodePairs(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,measure)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('measure','var');                 measure='phase';            end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPlots = getPlotHandles(2,1,[0.05 0.05 0.9 0.9],0.05,0.1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Check for intermediate data %%%%%%%%%%%%%%%%%%%
folderSavedData = fullfile(pwd,'savedData');
pairwiseDataToSave = fullfile(folderSavedData,['pairwiseMeanData' conditionType num2str(targetOnsetMatchingChoice) 'N' num2str(numTrialCutoff) 'TW' num2str(TWNum) measure '.mat']);

if exist(pairwiseDataToSave,'file')
    disp(['Loading saved data in ' pairwiseDataToSave]);
    load(pairwiseDataToSave,'pairwiseMeasureData','freqValsMT');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getMTValsSingleElectrode(TWNum,measure);
    
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
    %%%%%%%%%%%%%%%%%%%% Select sessions with numTrials>=cutoff %%%%%%%%%%%%%%%
    allNumTrialsMatrix = cell2mat(allNumTrials');
    minTrialsConditions = min(allNumTrialsMatrix(:,conditionsToUse),[],2);
    
    badSessionList = find(minTrialsConditions<=numTrialCutoff);
    goodSessionList = setdiff(1:length(allNumTrials),badSessionList);
    numGoodSessions = length(goodSessionList);
    disp(['Discarded sessions: ' num2str(badSessionList')]);
    
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
                        tmpMeasureList = cat(2,tmpMeasureList,getMeasure(d1,d2,measure));
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
numGoodSessions = size(pairwiseMeasureData,3);
for side=1:2
    for c=1:4 % Comparison 
        tmpIndexToCompare = indicesToCompare{c};
        disp(['side' num2str(side) '_' originalConditionsList{tmpIndexToCompare(1)} '-' originalConditionsList{tmpIndexToCompare(2)}]);
    
        meanDataAllSessions=[];
        for i=1:numGoodSessions
            tmpData1 = pairwiseMeasureData{side,tmpIndexToCompare(1),i};
            tmpData2 = pairwiseMeasureData{side,tmpIndexToCompare(2),i};
            meanDataAllSessions = cat(2,meanDataAllSessions,tmpData1-tmpData2);
        end
        allMeans0{side,c} = meanDataAllSessions;
    end
end

allMeans = combineComparisonDataAcrossBothArrays(allMeans0);

% Plot Means
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    plotData(hPlots(plotPos,1),freqValsMT,allMeans{i}',colorsForComparison{i},0);
end
% Plot SEMs
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    plotData(hPlots(plotPos,1),freqValsMT,allMeans{i}',colorsForComparison{i},1);
end
plot(hPlots(1,1),freqValsMT,zeros(1,length(freqValsMT)),'k--');
plot(hPlots(2,1),freqValsMT,zeros(1,length(freqValsMT)),'k--');
if strcmp(measure,'phase')
    ylabel(hPlots(1,1),'\Delta PPC');
else
    ylabel(hPlots(1,1),'\Delta Corr');
end
legend(hPlots(1,1),legendForComparison(1:2));
legend(hPlots(2,1),legendForComparison(3:4));
end
function outMeasure = getMeasure(d1,d2,measure)

numFreqs = size(d1,1);
outMeasure = zeros(numFreqs,1);
for i=1:numFreqs
    tmpData1 = squeeze(d1(i,:,:));
    tmpData2 = squeeze(d2(i,:,:));
    if strcmp(measure,'phase') % Get PPC
        outMeasure(i) = getPPC(tmpData1(:)-tmpData2(:));
    elseif strcmp(measure,'power') % Get Correlation
        cc = corrcoef(mean(tmpData1,1),mean(tmpData2,1));
        outMeasure(i) = cc(1,2);
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