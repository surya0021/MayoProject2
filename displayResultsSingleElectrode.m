% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis

function displayResultsSingleElectrode(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,showAUCFlag)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('showAUCFlag','var');             showAUCFlag=0;              end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNames = 'brcm'; % Attend-In Hit, Attend-Out Hit, Attend-In Miss, Attend-Out Miss: This is the order in which data is plotted
legendStrList = [{'AIH'} {'AOH'} {'AIM'} {'AOM'}];
colorsForComparison{1} = [{'k'} {[0.5 0.5 0.5]}]; % For comparison between AttIn and AttOut
colorsForComparison{2} = [{'b'} {'r'}]; % For comparison between Hits vs Misses
legendForComparison{1} = 'AI vs AO: H(black) & M(gray)';
legendForComparison{2} = 'H vs M: AI(b) & AO(r)';

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
if showAUCFlag
    numRows=5;
else
    numRows=4;
end
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
    plot(hPlots(1,1),c,h,colorNames(i)); hold(hPlots(1,1),'on');
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

for i=1:numGoodSessions
    allMeanFiringRates = cat(1,allMeanFiringRates,combineDataAcrossBothArrays(getMean(goodFiringRates{i}(:,conditionsToUse))));
    allMeanMTPower = cat(1,allMeanMTPower,combineDataAcrossBothArrays(getMean(goodMTPower{i}(:,conditionsToUse))));
end

% Firing Rates
makeBarPlot(hPlots(2,1),squeeze(allMeanFiringRates),colorNames,legendStrList);
title(hPlots(2,1),'Firing Rates'); ylabel(hPlots(2,1),'Spikes/s');

% MT Power
tmpData = 10*log10(allMeanMTPower); tmpFreq = freqValsMT; conditionStr='A';

for i=1:4
    plotData(hPlots(1,2),tmpFreq,squeeze(tmpData(:,:,i)),colorNames(i),conditionStr); % Change from AOV condition
    plotData(hPlots(2,2),tmpFreq,squeeze(tmpData(:,:,i) - tmpData(:,:,2)),colorNames(i),conditionStr); % Change from AOV condition
end

% Attention Difference
plotData(hPlots(1,3),tmpFreq,squeeze(tmpData(:,:,1) - tmpData(:,:,2)),colorsForComparison{1}{1},conditionStr);
plotData(hPlots(1,3),tmpFreq,squeeze(tmpData(:,:,3) - tmpData(:,:,4)),colorsForComparison{1}{2},conditionStr);
plot(hPlots(1,3),tmpFreq,zeros(1,length(tmpFreq)),'r');
title(hPlots(1,3),legendForComparison{1});

% Behavioral Difference
plotData(hPlots(2,3),tmpFreq,squeeze(tmpData(:,:,1) - tmpData(:,:,3)),colorsForComparison{2}{1},conditionStr);
plotData(hPlots(2,3),tmpFreq,squeeze(tmpData(:,:,2) - tmpData(:,:,4)),colorsForComparison{2}{2},conditionStr);
plot(hPlots(2,3),tmpFreq,zeros(1,length(tmpFreq)),'k');
title(hPlots(2,3),legendForComparison{2});

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot AUC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showAUCFlag
    numMeasuresToShow=2;
else
    numMeasuresToShow=1;
end
for measureType = 1:numMeasuresToShow % 1 - DPrime, 2 - AUC
    
    if measureType==1
        measureName = 'DPrime';
    elseif measureType==2
        measureName = 'AUC';
    end
    for comparisonType =  1:2    % 1 - Att In-out, 2 - H vs M
        allMeasuresFiringRates = [];
        allMeasuresMTPower = [];
        
        for i=1:numGoodSessions
            disp(['Getting ' measureName ' for session: ' num2str(i)]);
            allMeasuresFiringRates = cat(1,allMeasuresFiringRates,combineDataAcrossBothArrays(getAUC(goodFiringRates{i}(:,conditionsToUse),'A',comparisonType,measureType)));
            allMeasuresMTPower = cat(1,allMeasuresMTPower,combineDataAcrossBothArrays(getAUC(goodMTPower{i}(:,conditionsToUse),'A',comparisonType,measureType)));
        end
        
        % MT Power
        plotData(hPlots(comparisonType,3+measureType),freqValsMT,squeeze(allMeasuresMTPower(:,:,1)),colorsForComparison{comparisonType}{1});
        plotData(hPlots(comparisonType,3+measureType),freqValsMT,squeeze(allMeasuresMTPower(:,:,2)),colorsForComparison{comparisonType}{2});
        
        % Compare with AUC of firing rate
        plotData(hPlots(comparisonType,3+measureType),freqValsMT,repmat(squeeze(allMeasuresFiringRates(:,:,1)),1,length(freqValsMT)),colorsForComparison{comparisonType}{1});
        plotData(hPlots(comparisonType,3+measureType),freqValsMT,repmat(squeeze(allMeasuresFiringRates(:,:,2)),1,length(freqValsMT)),colorsForComparison{comparisonType}{2});
        
        title(hPlots(comparisonType,3+measureType),legendForComparison{comparisonType});
    end
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
function y=getAUC(x,condition,comparisonType,measureType)

if ~exist('condition','var');       condition = 'A';                    end
if ~exist('comparisonType','var');  comparisonType = 1;                 end
if ~exist('measureType','var');     measureType=1;                       end

if (size(x,1)~=2) || (size(x,2)~=4)
    error('Data matrix has inconsistent size');
end
xSize = numel(size(x{1,1})); % 2 for FR, 3 for LFP
y = cell(2,2);

for i=1:2
    for j=1:2
        if strcmp(condition,'A')
            if comparisonType==1
                d1 = abs(x{i,2*(j-1)+1}); d2 = abs(x{i,2*j}); % In - Out
            else
                d1 = abs(x{i,j}); d2 = abs(x{i,j+2}); % H vs M
            end
        elseif strcmp(condition,'P')
            if comparisonType==1
                d1 = angle(x{i,2*(j-1)+1}); d2 = angle(x{i,2*j}); % In - Out
            else
                d1 = angle(x{i,j}); d2 = angle(x{i,j+2}); % H vs M
            end
        end
        numElecs = size(d1,1);
        if size(d2,1)~=numElecs
            error('Number of electrodes are different');
        end
        
        if xSize==2 % Firing Rate
            tmpData = zeros(numElecs,1);
            for k=1:numElecs
                if measureType==1
                    tmpData(k) = getDPrime(d1(k,:),d2(k,:));               
                elseif measureType==2
                    tmpData(k) = ROCAnalysis(d1(k,:),d2(k,:));
                end
            end
        else % LFP power or phase
            numFreqs = size(d1,2);
            tmpData = zeros(numElecs,numFreqs);
            for k=1:numElecs
                for f=1:numFreqs
                    if measureType==1
                        tmpData(k,f) = getDPrime(squeeze(d1(k,f,:)),squeeze(d2(k,f,:)));
                    elseif measureType==2
                        tmpData(k,f) = ROCAnalysis(squeeze(d1(k,f,:)),squeeze(d2(k,f,:)));
                    end
                end
            end
        end
        y{i,j} = tmpData;
    end
end
end
function combinedData = combineDataAcrossBothArrays(data)

if (size(data,1)==2) && (size(data,2)==4)
    tmpCombinedData{1} = cat(1,data{1,1},data{2,2}); % Attend In Hit [(R)H0 & (L)H1]
    tmpCombinedData{2} = cat(1,data{1,2},data{2,1}); % Attend Out Hit [(R)H1 & (L)H0)]
    tmpCombinedData{3} = cat(1,data{1,3},data{2,4}); % Attend In Miss [(R)M0 & (L)M1]
    tmpCombinedData{4} = cat(1,data{1,4},data{2,3}); % Attend Out Miss [(R)M1 & (L)M0)
    
    combinedData = cat(3,tmpCombinedData{1},tmpCombinedData{2},tmpCombinedData{3},tmpCombinedData{4});
    
elseif (size(data,1)==2) && (size(data,2)==2) % For combining AUC results
    tmpCombinedData{1} = cat(1,data{1,1},data{2,2}); % Attend In AUC [(R)H0 & (L)H1]
    tmpCombinedData{2} = cat(1,data{1,2},data{2,1}); % Attend Out AUC [(R)H1 & (L)H0)]
    
    combinedData = cat(3,tmpCombinedData{1},tmpCombinedData{2});
end
end
function makeBarPlot(h,data,colorNames,legendStr)

N = size(data,1);

mData = mean(data,1);
semData = std(data,[],1)/sqrt(N);

for i=1:size(data,2)
    plot(h,i,mData(i),'color',colorNames(i),'marker','o');
    hold(h,'on');
    errorbar(h,i,mData(i),semData(i),'color',colorNames(i));
end
set(h,'XTick',1:4,'XTicklabel',legendStr);
xlim(h,[0 5]);
end
function plotData(hPlot,xs,data,colorName,condition)

if ~exist('condition','var');       condition = 'A';                    end

colorName2 = [0.5 0.5 0.5];

if strcmp(condition,'A')
    mData = squeeze(mean(data,1));
    sData = std(data,[],1)/sqrt(size(data,1));
    
    xsLong = [xs fliplr(xs)];
    ysLong = [mData+sData fliplr(mData-sData)];
    patch(xsLong,ysLong,colorName2,'parent',hPlot);
else
    mData = squeeze(circ_mean(data,[],1));
end
    
hold(hPlot,'on');
plot(hPlot,xs,mData,'color',colorName,'linewidth',1); 
end
function d = getDPrime(x1,x2)
stdVal = sqrt((var(x1)+var(x2))/2);
d = (mean(x1)- mean(x2))/stdVal;
end