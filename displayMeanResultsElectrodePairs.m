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
indicesToCompare{3} = [1 3]; % H vs M for Attend L
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
    load(pairwiseDataToSave,'pairwiseMeasureDataFR','pairwiseMeasureDataMT','freqValsMT');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [allFiringRates0,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getMTValsSingleElectrode(TWNum,measure);
    
    numSessions = length(allTargetOnsetTimes0);
    numConditions = length(allTargetOnsetTimes0{1});
    
    %%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    targetTimeBinWidthMS = 250;
    goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);
    
    %%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
    allFiringRates = cell(1,numSessions);
    allMTMeasure = cell(1,numSessions);
    allNumTrials = cell(1,numSessions);
    allTargetOnsetTimes = cell(1,numSessions);
    
    for i=1:numSessions
        tmpFiringRates = cell(2,numConditions);
        tmpMTMeasure = cell(2,numConditions);
        tmpAllNumTrials = zeros(1,numConditions);
        tmpTargetOnsetTimes = cell(1,numConditions);
        
        for k=1:numConditions
            for j=1:2
                tmpFiringRates{j,k} = allFiringRates0{i}{j,k}(:,:,:,goodStimNums{i}{k});
                tmpMTMeasure{j,k} = allMTMeasure0{i}{j,k}(:,:,:,goodStimNums{i}{k});
            end
            tmpAllNumTrials(k) = length(goodStimNums{i}{k});
            tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
        end
        allFiringRates{i} = tmpFiringRates;
        allMTMeasure{i} = tmpMTMeasure;
        allNumTrials{i} = tmpAllNumTrials;
        allTargetOnsetTimes{i} = tmpTargetOnsetTimes;
    end
    
    clear allFiringRates0 allMTMeasure0
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
    goodFiringRates = allFiringRates(goodSessionList);
    goodMTMeasure = allMTMeasure(goodSessionList);
    pairwiseMeasureDataFR = cell(2,numConditionsToUse,numGoodSessions);
    pairwiseMeasureDataMT = cell(2,numConditionsToUse,numGoodSessions);
    for side=1:2
        for c=1:numConditionsToUse
            for i=1:numGoodSessions
                disp(['side' num2str(side) '_' originalConditionsList{c} '_session' num2str(i)]);
                tmpDataFR0 = goodFiringRates{i}{side,conditionsToUse(c)};
                tmpDataMT = goodMTMeasure{i}{side,conditionsToUse(c)};
                [numElecs,numFreqs,numTapers,~] = size(tmpDataMT);
                % matching the dimensions of MT and firing rate (elec x
                % freq x tapers x trials)
                tmpDataFR = repmat(tmpDataFR0,[1 numFreqs numTapers 1]);
                
                tmpMeasureListFR = [];
                tmpMeasureListMT = [];
                for e1=1:numElecs
                    for e2=e1+1:numElecs
                        d1MT = squeeze(tmpDataMT(e1,:,:,:));
                        d2MT = squeeze(tmpDataMT(e2,:,:,:));
                        tmpMeasureListMT = cat(2,tmpMeasureListMT,getMeasure(d1MT,d2MT,measure));
                        
                        d1FR = squeeze(tmpDataFR(e1,:,:,:));
                        d2FR = squeeze(tmpDataFR(e2,:,:,:));
                        tmpMeasureListFR = cat(2,tmpMeasureListFR,getMeasure(d1FR,d2FR,'power'));
                    end
                end
                pairwiseMeasureDataFR{side,c,i} = tmpMeasureListFR;
                pairwiseMeasureDataMT{side,c,i} = tmpMeasureListMT;
            end
        end
    end
    save(pairwiseDataToSave,'pairwiseMeasureDataFR','pairwiseMeasureDataMT','freqValsMT');
end

%%%%%%%%%%%%%%%%%%%% Get Means and DPrimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pairwise data is for conditions H0, H1, M0 and M1. This needs to be
% converted to Hit(InVsOut), Miss(InVsOut), In(HitVsMiss) and Out(HitVsMiss)
for dataType=1:2  % MT measure and Firing rates
    clear allMeans pairwiseMeasureData
    if dataType==1
        pairwiseMeasureData = pairwiseMeasureDataMT;
    elseif dataType==2
        pairwiseMeasureData = pairwiseMeasureDataFR;
    end
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
    if dataType==1
        allMeansMT = allMeans;
    elseif dataType==2
        allMeansFR = allMeans;
    end
end
% Plot Means
pMT = zeros(1,4); pFR = zeros(1,4);
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    pMT(i) = plotData(hPlots(plotPos,1),freqValsMT,allMeansMT{i}',colorsForComparison{i},0);
    hold(hPlots(plotPos,1),'on')
    if strcmp(measure,'power') % Plot firing rate correlation with power correlation only
        pFR(i) = plotData(hPlots(plotPos,1),freqValsMT,allMeansFR{i}',colorsForComparison{i},1);
    end
end

plot(hPlots(1,1),freqValsMT,zeros(1,length(freqValsMT)),'k--');
plot(hPlots(2,1),freqValsMT,zeros(1,length(freqValsMT)),'k--');
if strcmp(measure,'phase')
    ylabel(hPlots(1,1),'\Delta PPC');
    xlabel(hPlots(1,1),'Frequency');
else
    ylabel(hPlots(1,1),'\Delta Corr');
    xlabel(hPlots(1,1),'Frequency');
end
if strcmp(measure,'power')
    legendForComparisonMT = cellfun(@(x) [x ' Power'],legendForComparison,'un',0);
    legendForComparisonFR = cellfun(@(x) [x ' Firing rate'],legendForComparison,'un',0);
    legend(hPlots(1,1),[pMT(1:2) pFR(1:2)],[legendForComparisonMT(1:2) legendForComparisonFR(1:2)]);
    legend(hPlots(2,1),[pMT(3:4) pFR(3:4)],[legendForComparisonMT(3:4) legendForComparisonFR(3:4)]);
else
    legend(hPlots(1,1),pMT(1:2),legendForComparison(1:2));
    legend(hPlots(2,1),pMT(3:4),legendForComparison(3:4));
end
end

function outMeasure = getMeasure(d1,d2,measure)

numFreqs = size(d1,1);
outMeasure = zeros(numFreqs,1);
for i=1:numFreqs
    tmpData1 = squeeze(d1(i,:,:));
    tmpData2 = squeeze(d2(i,:,:));
    if strcmp(measure,'phase') % Get PPC
        outMeasure(i) = getPPC(circ_mean(tmpData1-tmpData2,[],1));
    elseif strcmp(measure,'power') % Get Correlation
        cc = corrcoef(mean(tmpData1,1),mean(tmpData2,1));
        outMeasure(i) = cc(1,2);
    end
end
end
function combinedData = combineComparisonDataAcrossBothArrays(data)

combinedData = cell(1,4);

% Input data has size of 2x4. Rows correspond to array R and L respectively.
% The 4 comparisons are 1-2, 3-4, 1-3 and 3-4, corresponding to H0-H1,
% M0-M1, H0-M0 and H1-M1. Here 0 and 1 correspond to sides L and R.

% AttIn-AttOut Hit & Miss
% The comparison 1-2 (H0-H1) corresponds to AttIn-out Hit for array 1(R) and
% AttOut - AttIn Hit for array 2(L). Therefore We can get AttIn-Out Hit by
% simply flipping the values for array 2 and concatenating

combinedData{1} = [data{1,1} -data{2,1}];
combinedData{2} = [data{1,2} -data{2,2}];

% Hits-Miss for AttIn and AttOut: 
% For comparison 3 (1-3 or H0-M0) - these correspond to AttInHvsM for array
% 0 but AttOutHVsM for array 1. This flips for comparison 4.

combinedData{3} = [data{1,3} data{2,4}];
combinedData{4} = [data{1,4} data{2,3}];

end
function p = plotData(hPlot,xs,data,colorName,plotFRFlag)

if ~exist('plotFRFlag','var');     plotFRFlag = 1;                    end

if plotFRFlag
    lineType = '--';
else
    lineType = '-';
end
tmp=rgb2hsv(colorName); 
tmp2 = [tmp(1) tmp(2)/3 tmp(3)];
colorName2 = hsv2rgb(tmp2); % Same color with more saturation

mData = squeeze(mean(data,1));
sData = std(data,[],1)/sqrt(size(data,1));
xsLong = [xs fliplr(xs)];
ysLong = [mData+sData fliplr(mData-sData)];
patch(xsLong,ysLong,colorName2,'EdgeColor','none','parent',hPlot);

hold(hPlot,'on');
p = plot(hPlot,xs,mData,'color',colorName,'linewidth',2,'lineStyle',lineType); 
end