% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis

function [allMeanMTPPC,freqValsMT] = displayResultsMTPhase(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,displayFlag)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('displayFlag','var');             displayFlag=1;              end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNamesList(1,:) = [0 0 1]; % Attend-In Hit Blue
colorNamesList(2,:) = [1 0 0]; % Attend-Out Hit Red
colorNamesList(3,:) = [0 1 1]; % Attend-In Miss Cyan
colorNamesList(4,:) = [1 0 1]; % Attend-Out Miss: Magenta
legendStrList = [{'AIH'} {'AOH'} {'AIM'} {'AOM'}];

% To plot differences
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
[allFiringRates0,allMTPhase0,~,allTargetOnsetTimes0,freqValsMT] = getMTValsSingleElectrode(TWNum,'phase');
numSessions = length(allTargetOnsetTimes0);
numConditions = length(allTargetOnsetTimes0{1});

%%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetTimeBinWidthMS = 250;
goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);

%%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
allFiringRates = cell(1,numSessions);
allMTPhase = cell(1,numSessions);
allNumTrials = cell(1,numSessions);
allTargetOnsetTimes = cell(1,numSessions);

for i=1:numSessions
    tmpFiringRates = cell(2,numConditions);
    tmpMTPhase = cell(2,numConditions);
    tmpAllNumTrials = zeros(1,numConditions);
    tmpTargetOnsetTimes = cell(1,numConditions);
    
    for k=1:numConditions
        for j=1:2
            tmpFiringRates{j,k} = allFiringRates0{i}{j,k}(:,goodStimNums{i}{k});
            tmpMTPhase{j,k} = allMTPhase0{i}{j,k}(:,:,:,goodStimNums{i}{k});
        end
        tmpAllNumTrials(k) = length(goodStimNums{i}{k});
        tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
    end
    allFiringRates{i} = tmpFiringRates;
    allMTPhase{i} = tmpMTPhase;
    allNumTrials{i} = tmpAllNumTrials;
    allTargetOnsetTimes{i} = tmpTargetOnsetTimes;
end

clear allMTPhase0;
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
if displayFlag
    numRows=3;
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Mean Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

goodFiringRates = allFiringRates(goodSessionList);
goodMTPhase = allMTPhase(goodSessionList);

allMeanFiringRates = [];
allMeanMTPPC = [];

for i=1:numGoodSessions
    allMeanFiringRates = cat(1,allMeanFiringRates,combineDataAcrossBothArrays(getMean(goodFiringRates{i}(:,conditionsToUse))));
    allMeanMTPPC = cat(1,allMeanMTPPC,combineDataAcrossBothArrays(getMeanPPC(goodMTPhase{i}(:,conditionsToUse))));
end
if displayFlag
    % Firing Rates
    makeBarPlot(hPlots(2,1),squeeze(allMeanFiringRates),colorNamesList,legendStrList);
    title(hPlots(2,1),'Firing Rates'); ylabel(hPlots(2,1),'Spikes/s');
    
    % MT PPC
    tmpData = allMeanMTPPC; tmpFreq = freqValsMT;
    for i=1:4
        plotData(hPlots(1,2),tmpFreq,squeeze(tmpData(:,:,i)),colorNamesList(i,:));
        plotData(hPlots(2,2),tmpFreq,squeeze(tmpData(:,:,i) - tmpData(:,:,2)),colorNamesList(i,:)); % Change from AOV condition
    end
    ylabel(hPlots(1,2),'PPC'); plot(hPlots(1,2),tmpFreq,zeros(1,length(tmpFreq)),'k--');
    ylabel(hPlots(2,2),'\DeltaPPC');
    
    % Plot the differences
    % First plot without errobars
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        tmpIndexToCompare = indicesToCompare{i};
        plotData(hPlots(plotPos,3),tmpFreq,squeeze(tmpData(:,:,tmpIndexToCompare(1)) - tmpData(:,:,tmpIndexToCompare(2))),colorsForComparison{i},0);
    end
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        tmpIndexToCompare = indicesToCompare{i};
        plotData(hPlots(plotPos,3),tmpFreq,squeeze(tmpData(:,:,tmpIndexToCompare(1)) - tmpData(:,:,tmpIndexToCompare(2))),colorsForComparison{i});
    end
    
    plot(hPlots(1,3),tmpFreq,zeros(1,length(tmpFreq)),'k--');
    legend(hPlots(1,3),[legendForComparison(1) legendForComparison(2)],'location','best');
    ylabel(hPlots(1,3),'\DeltaPPC');
    plot(hPlots(2,3),tmpFreq,zeros(1,length(tmpFreq)),'k--');
    legend(hPlots(2,3),[legendForComparison(3) legendForComparison(4)],'location','best');
    ylabel(hPlots(2,3),'\DeltaPPC');
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
function y=getMeanPPC(x)

num1 = size(x,1); 
num2 = size(x,2);
y = cell(num1,num2);

for i=1:num1
    for j=1:num2
        allPhases = x{i,j};
        [numElectrodes,numFreqs,numTapers,~] = size(allPhases);
        
        mTMPPPC = zeros(numElectrodes,numFreqs);
        for e=1:numElectrodes
            for f=1:numFreqs
                tmpPPC = zeros(1,numTapers);
                for t=1:numTapers
                    tmpPPC(t) = getPPC(squeeze(allPhases(e,f,t,:)));
                end
                mTMPPPC(e,f) = mean(tmpPPC);
            end
        end
        y{i,j} = mTMPPPC;
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
end
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