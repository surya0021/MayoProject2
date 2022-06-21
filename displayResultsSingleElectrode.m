% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis

function displayResultsSingleElectrode(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNamesList(1,:) = [0 0 1]; % Attend-In Hit Blue
colorNamesList(2,:) = [1 0 0]; % Attend-Out Hit Red
colorNamesList(3,:) = [0 1 1]; % Attend-In Miss Cyan
colorNamesList(4,:) = [1 0 1]; % Attend-Out Miss: Magenta
legendStrList = [{'AIH'} {'AOH'} {'AIM'} {'AOM'}];

% To plot differences and dPrimes
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[allFiringRates0,allMTPower0,~,allTargetOnsetTimes0,freqValsMT] = getAnalysisMeasuresSingleElectrode(TWNum);
numSessions = length(allTargetOnsetTimes0);
numConditions = length(allTargetOnsetTimes0{1});

%%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetTimeBinWidthMS = 250;
goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);

%%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
allFiringRates = cell(1,numSessions);
allMTPower = cell(1,numSessions);
allNumTrials = cell(1,numSessions);
allTargetOnsetTimes = cell(1,numSessions);

for i=1:numSessions
    tmpFiringRates = cell(2,numConditions);
    tmpMTPower = cell(2,numConditions);
    tmpAllNumTrials = zeros(1,numConditions);
    tmpTargetOnsetTimes = cell(1,numConditions);
    
    for k=1:numConditions
        for j=1:2
            tmpFiringRates{j,k} = allFiringRates0{i}{j,k}(:,goodStimNums{i}{k});
            tmpMTPower{j,k} = allMTPower0{i}{j,k}(:,:,goodStimNums{i}{k});
        end
        tmpAllNumTrials(k) = length(goodStimNums{i}{k});
        tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
    end
    allFiringRates{i} = tmpFiringRates;
    allMTPower{i} = tmpMTPower;
    allNumTrials{i} = tmpAllNumTrials;
    allTargetOnsetTimes{i} = tmpTargetOnsetTimes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Condition Indices %%%%%%%%%%%%%%%%%%%%%%%%%%
% The order of the 12 conditions is as follows: {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'};
if strcmpi(conditionType(1),'V')
    conditionsToUse = [1 2 5 6];
elseif strcmpi(conditionType(1),'N')
    conditionsToUse = 9:12;
elseif strcmpi(conditionType(1),'I')
    conditionsToUse = [3 4 7 8];
end
numConditionsToUse = length(conditionsToUse);

%%%%%%%%%%%%%%%%%%%% Select sessions with numTrials>=cutoff %%%%%%%%%%%%%%%
allNumTrialsMatrix = cell2mat(allNumTrials');
minTrialsConditions = min(allNumTrialsMatrix(:,conditionsToUse),[],2);

badSessionList = find(minTrialsConditions<=numTrialCutoff);
goodSessionList = setdiff(1:length(allNumTrials),badSessionList);
numGoodSessions = length(goodSessionList);
disp(['Discarded sessions: ' num2str(badSessionList')]);

%%%%%%%%%%%%%%%%%%%%%%%%% Display Mean Responses %%%%%%%%%%%%%%%%%%%%%%%%%%
numRows=4;
hPlots = getPlotHandles(2,numRows,[0.05 0.05 0.9 0.9],0.05,0.1,0);

%%%%%%%%%%%%%%%%%%%%%%% Display TargetOnset Histogram %%%%%%%%%%%%%%%%%%%%%
targetOnsetTimesForHistogram = cell(1,numConditionsToUse);
for i=1:numGoodSessions
    for j=1:numConditionsToUse
        targetOnsetTimesForHistogram{j} = cat(2,targetOnsetTimesForHistogram{j},allTargetOnsetTimes{goodSessionList(i)}{conditionsToUse(j)});
    end
end

targetOnsetEdges = 500:targetTimeBinWidthMS:5500;
c = 500+targetTimeBinWidthMS/2:targetTimeBinWidthMS:5500;
legendStr = cell(1,numConditionsToUse);
for i=1:numConditionsToUse
    h = histcounts(targetOnsetTimesForHistogram{i},targetOnsetEdges);
    plot(hPlots(1,1),c,h,'color',colorNamesList(i,:)); hold(hPlots(1,1),'on');
    legendStr{i} = [legendStrList{i} '(' num2str(length(targetOnsetTimesForHistogram{i})) ')'];
end
legend(hPlots(1,1),legendStr);
xlabel(hPlots(1,1),'TargetOnset (ms)'); ylabel(hPlots(1,1),'Num Stim');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Mean Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodFiringRates = allFiringRates(goodSessionList);
goodMTPower = allMTPower(goodSessionList);

allMeanFiringRates = [];
allMeanMTPower = [];

% Input data is a 2x4 matrix (2 sides x 4 conditions - H0, H1, M0 and M1).
% They are combined across arrays to yield AttInHit, AttOutHit, AttInMiss and AttOutMiss
for i=1:numGoodSessions
    allMeanFiringRates = cat(1,allMeanFiringRates,combineDataAcrossBothArrays(getMean(goodFiringRates{i}(:,conditionsToUse))));
    allMeanMTPower = cat(1,allMeanMTPower,combineDataAcrossBothArrays(getMean(goodMTPower{i}(:,conditionsToUse))));
end

% Firing Rates
makeBarPlot(hPlots(2,1),squeeze(allMeanFiringRates),colorNamesList,legendStrList);
title(hPlots(2,1),'Firing Rates'); ylabel(hPlots(2,1),'Spikes/s');

% MT Power
tmpData = 10*log10(allMeanMTPower); tmpFreq = freqValsMT;

% Plot without SEM first to get the legend in proper order
for i=1:4
    plotData(hPlots(1,2),tmpFreq,squeeze(tmpData(:,:,i))/10,colorNamesList(i,:),0);
    plotData(hPlots(2,2),tmpFreq,squeeze(tmpData(:,:,i) - tmpData(:,:,2)),colorNamesList(i,:),0); % Change from AOV condition
end
% Now plot with SEMs
for i=1:4
    plotData(hPlots(1,2),tmpFreq,squeeze(tmpData(:,:,i))/10,colorNamesList(i,:),1);
    plotData(hPlots(2,2),tmpFreq,squeeze(tmpData(:,:,i) - tmpData(:,:,2)),colorNamesList(i,:),1); % Change from AOV condition
end
legend(hPlots(1,2),legendStrList); ylabel(hPlots(1,2),'Power');
legend(hPlots(2,2),legendStrList); ylabel(hPlots(2,2),'\DeltaPower (dB)');

% Plot the differences
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    tmpIndexToCompare = indicesToCompare{i};
    plotData(hPlots(plotPos,3),tmpFreq,squeeze(tmpData(:,:,tmpIndexToCompare(1)) - tmpData(:,:,tmpIndexToCompare(2))),colorsForComparison{i},0);
end
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    tmpIndexToCompare = indicesToCompare{i};
    plotData(hPlots(plotPos,3),tmpFreq,squeeze(tmpData(:,:,tmpIndexToCompare(1)) - tmpData(:,:,tmpIndexToCompare(2))),colorsForComparison{i},1);
end

plot(hPlots(1,3),tmpFreq,zeros(1,length(tmpFreq)),'k--');
legend(hPlots(1,3),[legendForComparison(1) legendForComparison(2)],'location','best');
ylabel(hPlots(1,3),'\DeltaPower (dB)');

plot(hPlots(2,3),tmpFreq,zeros(1,length(tmpFreq)),'k--');
legend(hPlots(2,3),[legendForComparison(3) legendForComparison(4)],'location','best');
ylabel(hPlots(2,3),'\DeltaPower (dB)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get DPrimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allDPrimesFiringRates0 = cell(4,numGoodSessions);
allDPrimesMTPower0 = cell(4,numGoodSessions);

for c=1:4 % Comparison
    tmpIndexToCompare = indicesToCompare{c};
    for i=1:numGoodSessions
        allDPrimesFiringRates0{c,i} = getDPrimeMatrix(goodFiringRates{i}(:,conditionsToUse(tmpIndexToCompare)));
        allDPrimesMTPower0{c,i} = getDPrimeMatrix(goodMTPower{i}(:,conditionsToUse(tmpIndexToCompare)));
    end
end

allDPrimesFiringRates = combineDPrimesAcrossBothArrays(allDPrimesFiringRates0);
allDPrimesMTPower = combineDPrimesAcrossBothArrays(allDPrimesMTPower0);

% Plot DPrimes
for i=1:4
    if (i<3); plotPos = 1; else; plotPos=2; end
    plotData(hPlots(plotPos,4),freqValsMT,allDPrimesMTPower{i},colorsForComparison{i},1);
    plotData(hPlots(plotPos,4),freqValsMT,repmat(squeeze(allDPrimesFiringRates{i}),1,length(freqValsMT)),colorsForComparison{i},1);
end
for i=1:2
    plot(hPlots(i,4),freqValsMT,zeros(1,length(freqValsMT)),'k--');
    ylabel(hPlots(i,4),'dPrime');
end
end
function y=getMean(x,condition)

if ~exist('condition','var');       condition = 'A';                    end

num1 = size(x,1); 
num2 = size(x,2);

xSize = numel(size(x{1,1}));
y = cell(num1,num2);

for i=1:num1
    for j=1:num2
        if strcmp(condition,'A') % Amplitude
            y{i,j} = mean(abs(x{i,j}),xSize);
        elseif strcmp(condition,'P') % Phase
            y{i,j} = circ_mean(angle(x{i,j}),[],xSize);
        end
    end
end
end
function y=getDPrimeMatrix(x)

if (size(x,1)~=2) || (size(x,2)~=2)
    error('Data matrix has inconsistent size');
end
xSize = numel(size(x{1,1})); % 2 for FR, 3 for LFP

y = cell(1,2);
for i=1:2
    d1 = x{i,1}; d2 = x{i,2};
    
    numElecs = size(d1,1);
    if size(d2,1)~=numElecs
        error('Number of electrodes are different');
    end
    
    if xSize==2 % Firing Rate
        tmpData = zeros(numElecs,1);
        for k=1:numElecs
            tmpData(k) = getDPrime(d1(k,:),d2(k,:));
        end
    else % LFP power or phase
        numFreqs = size(d1,2);
        tmpData = zeros(numElecs,numFreqs);
        for k=1:numElecs
            for f=1:numFreqs
                tmpData(k,f) = getDPrime(squeeze(d1(k,f,:)),squeeze(d2(k,f,:)));
            end
        end
    end
    y{i} = tmpData;
end
end
function combinedData = combineDataAcrossBothArrays(data)

if (size(data,1)==2) && (size(data,2)==4)
    tmpCombinedData{1} = cat(1,data{1,1},data{2,2}); % Attend In Hit [(R)H0 & (L)H1]
    tmpCombinedData{2} = cat(1,data{1,2},data{2,1}); % Attend Out Hit [(R)H1 & (L)H0)]
    tmpCombinedData{3} = cat(1,data{1,3},data{2,4}); % Attend In Miss [(R)M0 & (L)M1]
    tmpCombinedData{4} = cat(1,data{1,4},data{2,3}); % Attend Out Miss [(R)M1 & (L)M0)
    
    combinedData = cat(3,tmpCombinedData{1},tmpCombinedData{2},tmpCombinedData{3},tmpCombinedData{4});
end
end
function combinedData = combineDPrimesAcrossBothArrays(data)

numGoodSessions = size(data,2);
combinedData = cell(1,4);

% Input data has size of 4 x numGoodSessions. Each entry is a cell array of
% size 2x1, each containing data from one array. Data needs to be combined
% across arrays as follows:

% The 4 comparisons are 1-2, 3-4, 1-3 and 3-4, corresponding to H0-H1,
% M0-M1, H0-M0 and H1-M1. Here 0 and 1 correspond to sides 0 and 1.

% AttIn-AttOut Hit & Miss
% The comparison 1-2 (H0-H1) corresponds to AttIn-out Hit for side 0 and
% AttOut - AttIn Hit for side 1. Therefore We can get AttIn-Out Hit by
% simply flipping the values for side 1 and concatenating

for c=1:2
    allVals = [];
    for i=1:numGoodSessions
        tmpData = data{c,i};
        allVals = cat(1,allVals,tmpData{1});
        allVals = cat(1,allVals,-tmpData{2});
    end
    combinedData{c} = allVals;
end

% Hits-Miss for AttIn and AttOut: 
% For comparison 3 (1-3 or H0-M0) - these correspond to AttInHvsM for array
% 0 but AttOutHVsM for array 1. This flips for comparison 4.

attInHitVsMiss = []; 
attOutHitVsMiss = [];
for i=1:numGoodSessions
    comparison3 = data{3,i};
    comparison4 = data{4,i};
    attInHitVsMiss = cat(1,attInHitVsMiss,comparison3{1});
    attInHitVsMiss = cat(1,attInHitVsMiss,comparison4{2});
    attOutHitVsMiss = cat(1,attOutHitVsMiss,comparison3{2});
    attOutHitVsMiss = cat(1,attOutHitVsMiss,comparison4{1});
end
combinedData{3} = attInHitVsMiss;
combinedData{4} = attOutHitVsMiss;

end
function makeBarPlot(h,data,colorNames,legendStr)

N = size(data,1);

mData = mean(data,1);
semData = std(data,[],1)/sqrt(N);

for i=1:size(data,2)
    plot(h,i,mData(i),'color',colorNames(i,:),'marker','o');
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color',colorNames(i,:));
end
set(h,'XTick',1:4,'XTicklabel',legendStr);
xlim(h,[0 5]);
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
function d = getDPrime(x1,x2)
stdVal = sqrt((var(x1)+var(x2))/2);
d = (mean(x1)- mean(x2))/stdVal;
end