% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis. If using FFT, this should be equal to numDivisions
% measure - phase or power
% useFFTFlag - uses FFT instead of MT

function displayResultsElectrodePairs(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,measure,useFFTFlag)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('measure','var');                 measure='phase';            end
if ~exist('useFFTFlag','var');              useFFTFlag=0;               end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the original conditions - AttLeftHit AttRightHit AttLeftMiss AttRightMiss
% colorNamesList(1,:) = [0 0 1]; % Attend-In Hit Blue
% colorNamesList(2,:) = [1 0 0]; % Attend-Out Hit Red
% colorNamesList(3,:) = [0 1 1]; % Attend-In Miss Cyan
% colorNamesList(4,:) = [1 0 1]; % Attend-Out Miss: Magenta

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
hPlots = getPlotHandles(2,2,[0.05 0.05 0.9 0.9],0.05,0.1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Check for intermediate data %%%%%%%%%%%%%%%%%%%
folderSavedData = fullfile(pwd,'savedData');

if useFFTFlag
    pairwiseDataToSave = fullfile(folderSavedData,['pairwiseDataFFT' conditionType num2str(targetOnsetMatchingChoice) 'N' num2str(numTrialCutoff) 'TW' num2str(TWNum) measure '.mat']);
else
    pairwiseDataToSave = fullfile(folderSavedData,['pairwiseData' conditionType num2str(targetOnsetMatchingChoice) 'N' num2str(numTrialCutoff) 'TW' num2str(TWNum) measure '.mat']);
end

if exist(pairwiseDataToSave,'file')
    disp(['Loading saved data in ' pairwiseDataToSave]);
    load(pairwiseDataToSave,'pairwiseMeasureDataMT','pairwiseMeasureDataFR','freqValsMT');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if useFFTFlag
        [allFiringRates0,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getFFTValsSingleElectrode(TWNum,measure); % Getting FFT Measures but calling them MT
    else
        [allFiringRates0,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getMTValsSingleElectrode(TWNum,measure);
    end
    
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
                tmpDataMT = goodMTMeasure{i}{side,conditionsToUse(c)};
                tmpDataFR0 = goodFiringRates{i}{side,conditionsToUse(c)};
                
                % making the dimensions of tmpData(MTVals) and tmpDataFR0(firing rates) equal.
                
                sizeMT = size(tmpDataMT); sizeFR = size(tmpDataFR0);
                repDim = [1 1 1 1];
                if sizeMT(2)~=sizeFR(2) && sizeMT(3)==sizeFR(3)  % when freqVals don't match but divisions/taper match.
                   repDim(2) = sizeMT(2); 
                elseif sizeMT(2)~=sizeFR(2) && sizeMT(3)~=sizeFR(3) % when freqVals and divisions/taper don't match.
                    repDim(2) = sizeMT(2);
                    repDim(3) = sizeMT(3);
                end
                tmpDataFR = repmat(tmpDataFR0,repDim);
                if ~isequal(sizeMT,size(tmpDataFR))
                    error('Dimensions of MT and firing rate matrix does not match')
                end
                
                numElecs = size(tmpDataMT,1);
                
                tmpMeasureListMT = [];
                tmpMeasureListFR = [];
                for e1=1:numElecs
                    for e2=e1+1:numElecs
                        d1MT = squeeze(tmpDataMT(e1,:,:,:));    
                        d2MT = squeeze(tmpDataMT(e2,:,:,:));
                        tmpMeasureListMT = cat(3,tmpMeasureListMT,getMeasure(d1MT,d2MT,measure));
                        if useFFTFlag==1
                            d1FR = squeeze(tmpDataFR(e1,:,:,:));
                            d2FR = squeeze(tmpDataFR(e2,:,:,:));
                            tmpMeasureListFR = cat(3,tmpMeasureListFR,getMeasure(d1FR,d2FR,'power'));
                        end
                    end
                end
                pairwiseMeasureDataMT{side,c,i} = tmpMeasureListMT;
                pairwiseMeasureDataFR{side,c,i} = tmpMeasureListFR;
            end
        end
    end
    save(pairwiseDataToSave,'pairwiseMeasureDataMT','pairwiseMeasureDataFR','freqValsMT');
end

%%%%%%%%%%%%%%%%%%%% Get Means and DPrimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pairwise data is for conditions H0, H1, M0 and M1. This needs to be
% converted to Hit(InVsOut), Miss(InVsOut), In(HitVsMiss) and Out(HitVsMiss)
if useFFTFlag==1
    nType=2;
else
    nType=1;
end
for dataType=1:nType
    clear allMeans allDPrimes
    if dataType==1
        pairwiseMeasureData = pairwiseMeasureDataMT;
    elseif dataType==2
        pairwiseMeasureData = pairwiseMeasureDataFR;
    end
    allMeans0 = cell(2,4);
    allDPrimes0 = cell(2,4);
    numGoodSessions = size(pairwiseMeasureData,3);
    for side=1:2
        for c=1:4 % Comparison
            tmpIndexToCompare = indicesToCompare{c};
            disp(['side' num2str(side) '_' originalConditionsList{tmpIndexToCompare(1)} '-' originalConditionsList{tmpIndexToCompare(2)}]);
            
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
                        tmpMeanData(f,p) = mean(x1,'omitnan')-mean(x2,'omitnan');
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
    if dataType==1
        allMeansMT = allMeans;
        allDPrimesMT = allDPrimes;
    elseif dataType==2
        allMeansFR = allMeans;
        allDPrimesFR = allDPrimes;
    end
end
% Plot Means
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    pMT(i,1) = plotData(hPlots(plotPos,1),freqValsMT,allMeansMT{i}',colorsForComparison{i},0); %#ok<*AGROW>
    pMT(i,2) = plotData(hPlots(plotPos,2),freqValsMT,allDPrimesMT{i}',colorsForComparison{i},0);
end

for i=1:2
    plot(hPlots(i,1),freqValsMT,zeros(1,length(freqValsMT)),'k--');
    plot(hPlots(i,2),freqValsMT,zeros(1,length(freqValsMT)),'k--');
    if strcmp(measure,'phase')
        ylabel(hPlots(i,1),'\Delta PPC');
    else
        ylabel(hPlots(i,1),'\Delta Corr');
    end
    ylabel(hPlots(i,2),'dPrime');
end
%%%%%%%%%%%%%%%%%%%% Plotting firing rate correlation %%%%%%%%%%%%%%%%%%%%%
pFR = [];
if useFFTFlag==1
    % Plot Means
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        if strcmp(measure,'power') %plot mean change in firing rate correlation with power only
            pFR(i,1) = plotData(hPlots(plotPos,1),freqValsMT,allMeansFR{i}',colorsForComparison{i},1);
        end
        pFR(i,2) = plotData(hPlots(plotPos,2),freqValsMT,allDPrimesFR{i}',colorsForComparison{i},1);
    end
    
end
if useFFTFlag==1 % Condition where firing rate correlations are plotted
    legendForComparisonMTdPrime = cellfun(@(x) [x ' ' measure],legendForComparison,'un',0);
    legendForComparisonFRdPrime = cellfun(@(x) [x ' Firing rate'],legendForComparison,'un',0);
    
    legend(hPlots(1,2),[pMT(1:2,2) ; pFR(1:2,2)],[legendForComparisonMTdPrime(1:2) legendForComparisonFRdPrime(1:2)]);
    legend(hPlots(2,2),[pMT(3:4,2) ; pFR(3:4,2)],[legendForComparisonMTdPrime(3:4) legendForComparisonFRdPrime(3:4)]);
    
    if strcmp(measure,'power')
        legend(hPlots(1,1),[pMT(1:2,1) ; pFR(1:2,1)],[legendForComparisonMTdPrime(1:2) legendForComparisonFRdPrime(1:2)]);
        legend(hPlots(2,1),[pMT(3:4,1) ; pFR(3:4,1)],[legendForComparisonMTdPrime(3:4) legendForComparisonFRdPrime(3:4)]);
    elseif strcmp(measure,'phase')
        legend(hPlots(1,1),pMT(1:2,1),legendForComparison(1:2));
        legend(hPlots(2,1),pMT(3:4,1),legendForComparison(3:4));
    end
else
    for i=1:2
        legend(hPlots(1,i),pMT(1:2,i),legendForComparison(1:2));
        legend(hPlots(2,i),pMT(3:4,i),legendForComparison(3:4));
    end
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
            cc = corrcoef(tmpData1,tmpData2);
            outMeasure(i,t) = cc(1,2);
        end
    end
end
end
function combinedData = combineComparisonDataAcrossBothArrays(data)

combinedData = cell(1,4);

% Input data has size of 2x4. Rows correspond to array R and L respectively.
% The 4 comparisons are 1-2, 3-4, 1-3 and 2-4, corresponding to H0-H1,
% M0-M1, H0-M0 and H1-M1. Here 0 and 1 correspond to sides L and R.

% AttIn-AttOut Hit & Miss
% The comparison 1-2 (H0-H1) corresponds to AttIn-out Hit for array 1 (R) and
% AttOut - AttIn Hit for array 2 (L). Therefore We can get AttIn-Out Hit by
% simply flipping the values for array 2(L) and concatenating

combinedData{1} = [data{1,1} -data{2,1}];
combinedData{2} = [data{1,2} -data{2,2}];

% Hits-Miss for AttIn and AttOut: 
% For comparison 3 (1-3 or H0-M0) - these correspond to AttInHvsM for array
% 1 (R) but AttOutHVsM for array 2 (L). This flips for comparison 4.

combinedData{3} = [data{1,3} data{2,4}];
combinedData{4} = [data{1,4} data{2,3}];

end
function p = plotData(hPlot,xs,data,colorName,plotFRFlag)

if ~exist('plotFRFlag','var');     plotFRFlag = 0;                    end
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
p = plot(hPlot,xs,mData,'color',colorName,'linestyle',lineType,'linewidth',2); 
end