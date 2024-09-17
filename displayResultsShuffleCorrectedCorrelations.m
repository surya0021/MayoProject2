% 05/12/23 SSP
% Adapted from displayResultsElectrodePairs.m to compute shuffle 
% corrections to the binwise correlation of firing rate and power. This
% step is incorporated in the nested function getMeasures
% 
% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - numDivisions
% measure - always set to 'power'


function [allMeansMT,allDPrimesMT,allMeansFR,allDPrimesFR,allAbsMeansMT,allAbsMeansFR] = displayResultsShuffleCorrectedCorrelations(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,measure,displayFlag)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=10;                    end
if ~exist('measure','var');                 measure='power';            end
if ~exist('displayFlag','var');             displayFlag=1;              end
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
if displayFlag
    hPlots = getPlotHandles(2,2,[0.05 0.05 0.9 0.9],0.05,0.1,0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% Check for intermediate data %%%%%%%%%%%%%%%%%%%
folderSavedData = fullfile(pwd,'savedData');
pairwiseDataToSave = fullfile(folderSavedData,['pairwiseDataCorrectedFFT' conditionType num2str(targetOnsetMatchingChoice) 'N' num2str(numTrialCutoff) 'TW' num2str(TWNum) measure '.mat']);

if exist(pairwiseDataToSave,'file')
    disp(['Loading saved data in ' pairwiseDataToSave]);
    load(pairwiseDataToSave,'pairwiseMeasureDataMT','pairwiseMeasureDataFR','freqValsMT');
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [allFiringRates0,allMTMeasure0,~,allTargetOnsetTimes0,freqValsMT] = getFFTValsSingleElectrode(TWNum,measure); % Getting FFT Measures but calling them MT
    
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
                
                
                % making the dimensions of tmpDataMT(FFTVals) and tmpDataFR0(firing rates) equal.
                sizeMT = size(tmpDataMT); sizeFR = size(tmpDataFR0);
                repDim = [1 1 1 1];
                if prod(sizeMT([1 3:4])==sizeFR([1 3:4]))  % Make sure all dimensions except 2 match
                    repDim(2) = sizeMT(2);
                else
                    error('Dimensions of MT and firing rate matrix does not match');
                end
                tmpDataFR = repmat(tmpDataFR0,repDim);
               
                
                numElecs = size(tmpDataMT,1);
                
                tmpMeasureListMT = [];
                tmpMeasureListFR = [];
                for e1=1:numElecs
                    for e2=e1+1:numElecs
                        d1MT = squeeze(tmpDataMT(e1,:,:,:));    
                        d2MT = squeeze(tmpDataMT(e2,:,:,:));
                        tmpMeasureListMT = cat(3,tmpMeasureListMT,getMeasure(d1MT,d2MT,measure));
                        
                        d1FR = squeeze(tmpDataFR(e1,:,:,:));
                        d2FR = squeeze(tmpDataFR(e2,:,:,:));
                        tmpMeasureListFR = cat(3,tmpMeasureListFR,getMeasure(d1FR,d2FR,measure));
                        
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
allMeansMT = [];    allDPrimesMT = [];      allAbsMeansMT = [];
allMeansFR = [];    allDPrimesFR = [];      allAbsMeansFR = [];

nType=2;

for dataType=1:nType
    clear allMeans allDPrimes allAbsMeans
    if dataType==1
        pairwiseMeasureData = pairwiseMeasureDataMT;
    elseif dataType==2
        pairwiseMeasureData = pairwiseMeasureDataFR;
    end
    allMeans0 = cell(2,4);
    allDPrimes0 = cell(2,4);
    allAbsMeans0 = cell(2,4);
    numGoodSessions = size(pairwiseMeasureData,3);
    for side=1:2
        for c=1:4 % Comparison
            tmpIndexToCompare = indicesToCompare{c};
            disp(['side' num2str(side) '_' originalConditionsList{tmpIndexToCompare(1)} '-' originalConditionsList{tmpIndexToCompare(2)}]);
            
            meanDataAllSessions=[];
            dPrimeDataAllSessions=[];
            absMeanDataAllSessions=[];
            for i=1:numGoodSessions
                tmpData0 = pairwiseMeasureData{side,c,i};
                tmpData1 = pairwiseMeasureData{side,tmpIndexToCompare(1),i};
                tmpData2 = pairwiseMeasureData{side,tmpIndexToCompare(2),i};
                
                [numFreqs,~,numPairs,numCorrTypes] = size(tmpData1);
                
                clear tmpMeanData tmpDPrime tmpAbsMeanData
                tmpMeanData = zeros(numFreqs,numPairs,numCorrTypes);
                tmpDPrimeData = zeros(numFreqs,numPairs,numCorrTypes);
                tmpAbsMeanData = zeros(numFreqs,numPairs,numCorrTypes);
                for f=1:numFreqs
                    for p=1:numPairs
                        x0 = squeeze(tmpData0(f,:,p,:));
                        x1 = squeeze(tmpData1(f,:,p,:));
                        x2 = squeeze(tmpData2(f,:,p,:));
                        tmpMeanData(f,p,:) = mean(x1,'omitnan')-mean(x2,'omitnan');
                        tmpDPrimeData(f,p,:) = getDPrime(x1,x2);
                        tmpAbsMeanData(f,p,:) = mean(x0,'omitnan');
                    end
                end
                meanDataAllSessions = cat(2,meanDataAllSessions,tmpMeanData);
                dPrimeDataAllSessions = cat(2,dPrimeDataAllSessions,tmpDPrimeData);
                absMeanDataAllSessions = cat(2,absMeanDataAllSessions,tmpAbsMeanData);
            end
            allMeans0{side,c} = meanDataAllSessions;
            allDPrimes0{side,c} = dPrimeDataAllSessions;
            allAbsMeans0{side,c} = absMeanDataAllSessions;
        end
    end
    
    allMeans = combineComparisonDataAcrossBothArrays(allMeans0);
    allDPrimes = combineComparisonDataAcrossBothArrays(allDPrimes0);
    allAbsMeans =  combineDataAcrossBothArrays(allAbsMeans0);
    if dataType==1
        allMeansMT = cat(4,allMeans{:});
        allDPrimesMT = cat(4,allDPrimes{:});
        allAbsMeansMT = cat(4,allAbsMeans{:});

        allMeansMT = permute(allMeansMT,[1 2 4 3]);
        allDPrimesMT = permute(allDPrimesMT,[1 2 4 3]);
        allAbsMeansMT = permute(allAbsMeansMT,[1 2 4 3]);
    elseif dataType==2
        allMeansFR = cat(4,allMeans{:});
        allDPrimesFR = cat(4,allDPrimes{:});
        allAbsMeansFR = cat(4,allAbsMeans{:});

        allMeansFR = permute(allMeansFR,[1 2 4 3]);
        allDPrimesFR = permute(allDPrimesFR,[1 2 4 3]);
        allAbsMeansFR = permute(allAbsMeansFR,[1 2 4 3]);
    end
end

if displayFlag
    % Plot Means
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        pMT(i,1) = plotData(hPlots(plotPos,1),freqValsMT,allMeansMT(:,:,i,1)',colorsForComparison{i},0); %#ok<*AGROW>
        pMT(i,2) = plotData(hPlots(plotPos,2),freqValsMT,allDPrimesMT(:,:,i,1)',colorsForComparison{i},0);
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
    
    %Plot Means
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        if strcmp(measure,'power') %plot mean change in firing rate correlation with power only
            pFR(i,1) = plotData(hPlots(plotPos,1),freqValsMT,allMeansFR(:,:,i,1)',colorsForComparison{i},1);
        end
        pFR(i,2) = plotData(hPlots(plotPos,2),freqValsMT,allDPrimesFR(:,:,i,1)',colorsForComparison{i},1);
    end

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
    
end
end
function rValues = getMeasure(d1,d2,measure)

[numFreqs,~,numTrials] = size(d1);

rValues = zeros(numFreqs,numTrials,1,3); % Last dimension (3) is for corrected,raw and shuffled correlation values respectively
for i=1:numFreqs
    tmpData1 = squeeze(d1(i,:,:));
    tmpData2 = squeeze(d2(i,:,:));
    if strcmp(measure,'power') % Get Correlation
        cc = corr(tmpData1,tmpData2);
        ind = 1:size(cc,1)+1:size(cc,1)^2; % indices for diagonal elements
        rValues(i,:,1,2) = diag(cc); % uncorrected(raw) correlation
        if isequaln(cc(ind),rValues(i,:,1,2))
            cc(ind)=nan;
        else
            error('Mismatch between diagonal elements and its index')
        end
        rValues(i,:,1,3) = mean(cc,2,'omitnan'); % shuffle correlation across non-simultaneous trials 
        rValues(i,:,1,1) = rValues(i,:,1,2)-rValues(i,:,1,3); % Corrected correlation
    end 
end
% setting corrected correlation values to -1 or 1 when they exceed the
% these limit.
rValues(rValues>1)=1;
rValues(rValues<-1)=-1;
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

function combinedData = combineDataAcrossBothArrays(data)

combinedData = cell(1,4);

% Input data has size of 2x4. Rows correspond to array R and L respectively.


% AttIn and AttOut for Hit


combinedData{1} = [data{1,1} data{2,2}]; % [HLV(R) HRV(L)] Attend-in Hit
combinedData{2} = [data{1,2} data{2,1}]; % [HRV(R) HLV(L)] Attend-out Hit

% AttIn and AttOut for Miss 

combinedData{3} = [data{1,3} data{2,4}]; % [MLV(R) MRV(L)] % Attend-in Miss
combinedData{4} = [data{1,4} data{2,3}]; % [MRV(R) MLV(L)] % Attend-out Miss

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