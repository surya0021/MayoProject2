% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)

function displayResultsSingleElectrode(conditionType,targetOnsetMatchingChoice,numTrialCutoff)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorNames = 'brcm'; % Attend-In Hit, Attend-Out Hit, Attend-In Miss, Attend-Out Miss: This is the order in which data is plotted
legendStrList = [{'AIH'} {'AOH'} {'AIM'} {'AOM'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[allFiringRates0,allMTPower0,allFFTVals0,~,allTargetOnsetTimes0,freqValsMT,freqValsFFT] = getAnalysisMeasuresSingleElectrode;
numSessions = length(allTargetOnsetTimes0);
numConditions = length(allTargetOnsetTimes0{1});

%%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetTimeBinWidthMS = 250;
goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);

%%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
allFiringRates = cell(1,numSessions);
allMTPower = cell(1,numSessions);
allFFTVals = cell(1,numSessions);
allNumTrials = cell(1,numSessions);
allTargetOnsetTimes = cell(1,numSessions);

for i=1:numSessions
    tmpFiringRates = cell(2,numConditions);
    tmpMTPower = cell(2,numConditions);
    tmpFFTVals = cell(2,numConditions);
    tmpAllNumTrials = zeros(1,numConditions);
    tmpTargetOnsetTimes = cell(1,numConditions);
    
    for k=1:numConditions
        for j=1:2
            tmpFiringRates{j,k} = allFiringRates0{i}{j,k}(:,goodStimNums{i}{k});
            tmpMTPower{j,k} = allMTPower0{i}{j,k}(:,:,goodStimNums{i}{k});
            tmpFFTVals{j,k} = allFFTVals0{i}{j,k}(:,:,goodStimNums{i}{k});
        end
        tmpAllNumTrials(k) = length(goodStimNums{i}{k});
        tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
    end
    allFiringRates{i} = tmpFiringRates;
    allMTPower{i} = tmpMTPower;
    allFFTVals{i} = tmpFFTVals;
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

if targetOnsetMatchingChoice==1
    hPlots = getPlotHandles(5,3,[0.05 0.05 0.4 0.9],0.05,0.05,0);
else
    hPlots = getPlotHandles(5,3,[0.55 0.05 0.4 0.9],0.05,0.05,0);
end

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
goodFFTVals = allFFTVals(goodSessionList);

allMeanFiringRates = [];
allMeanMTPower = [];
allMeanFFTAmplitude = [];
allMeanFFTPhase = [];

for i=1:numGoodSessions
    allMeanFiringRates = cat(1,allMeanFiringRates,combineDataAcrossBothArrays(getMean(goodFiringRates{i}(:,conditionsToUse))));
    allMeanMTPower = cat(1,allMeanMTPower,combineDataAcrossBothArrays(getMean(goodMTPower{i}(:,conditionsToUse))));
    allMeanFFTAmplitude = cat(1,allMeanFFTAmplitude,combineDataAcrossBothArrays(getMean(goodFFTVals{i}(:,conditionsToUse),'A')));
    allMeanFFTPhase = cat(1,allMeanFFTPhase,combineDataAcrossBothArrays(getMean(goodFFTVals{i}(:,conditionsToUse),'P')));
end

% Firing Rates
makeBarPlot(hPlots(1,2),squeeze(allMeanFiringRates),colorNames,legendStrList);
title(hPlots(1,2),'Firing Rates'); ylabel(hPlots(1,2),'Spikes/s');

% FFT/MT Power and phase
for pos=1:3
    if pos==1 % MTPower
        tmpData = 10*log10(allMeanMTPower); tmpFreq = freqValsMT; conditionStr='A';
    elseif pos==2 % FFT Amplitude
        tmpData = 20*log10(allMeanFFTAmplitude); tmpFreq = freqValsFFT(1:size(tmpData,2)); conditionStr='A';
    elseif pos==3 % FFt Phase
        tmpData = allMeanFFTPhase; tmpFreq = freqValsFFT(1:size(tmpData,2)); conditionStr='P';
    end
    
    % Change from AOV condition
    for i=1:4
        plotData(hPlots(2,pos),tmpFreq,squeeze(tmpData(:,:,i) - tmpData(:,:,2)),colorNames(i),conditionStr);
    end
    
    % Attention Difference
    plotData(hPlots(3,pos),tmpFreq,squeeze(tmpData(:,:,1) - tmpData(:,:,2)),'k',conditionStr);
    plotData(hPlots(3,pos),tmpFreq,squeeze(tmpData(:,:,3) - tmpData(:,:,4)),[0.5 0.5 0.5],conditionStr);
    plot(hPlots(3,pos),tmpFreq,zeros(1,length(tmpFreq)),'r');
    title(hPlots(3,pos),'AI vs AO: H(black) & M(gray)');
    
    % Behavioral Difference
    plotData(hPlots(4,pos),tmpFreq,squeeze(tmpData(:,:,1) - tmpData(:,:,3)),colorNames(1),conditionStr);
    plotData(hPlots(4,pos),tmpFreq,squeeze(tmpData(:,:,2) - tmpData(:,:,4)),colorNames(2),conditionStr);
    plot(hPlots(4,pos),tmpFreq,zeros(1,length(tmpFreq)),'k');
    title(hPlots(4,pos),'H vs M: AI(b) & AO(r)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot AUC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allAUCFiringRates = [];
allAUCMTPower = [];
allAUCFFTAmplitude = [];
allAUCFFTPhase = [];

for i=1:numGoodSessions
    disp(['Getting AUC/dPrime for session: ' num2str(i)]); 
    allAUCFiringRates = cat(1,allAUCFiringRates,combineDataAcrossBothArrays(getAUC(goodFiringRates{i}(:,conditionsToUse))));
    allAUCMTPower = cat(1,allAUCMTPower,combineDataAcrossBothArrays(getAUC(goodMTPower{i}(:,conditionsToUse))));
    allAUCFFTAmplitude = cat(1,allAUCFFTAmplitude,combineDataAcrossBothArrays(getAUC(goodFFTVals{i}(:,conditionsToUse),'A')));
    allAUCFFTPhase = cat(1,allAUCFFTPhase,combineDataAcrossBothArrays(getAUC(goodFFTVals{i}(:,conditionsToUse),'P')));
end

for pos=1:3
    if pos==1 % MTPower
        tmpData = allAUCMTPower; tmpFreq = freqValsMT;
    elseif pos==2 % FFT Amplitude
        tmpData = allAUCFFTAmplitude; tmpFreq = freqValsFFT(1:size(tmpData,2));
    elseif pos==3 % FFt Phase
        tmpData = allAUCFFTPhase; tmpFreq = freqValsFFT(1:size(tmpData,2));
    end
    
    % Behavioral Difference
    plotData(hPlots(5,pos),tmpFreq,squeeze(tmpData(:,:,1)),colorNames(1));
    plotData(hPlots(5,pos),tmpFreq,squeeze(tmpData(:,:,2)),colorNames(2));
    
    % Compare with AUC of firing rate
    plotData(hPlots(5,pos),tmpFreq,repmat(squeeze(allAUCFiringRates(:,:,1)),1,length(tmpFreq)),colorNames(1));
    plotData(hPlots(5,pos),tmpFreq,repmat(squeeze(allAUCFiringRates(:,:,2)),1,length(tmpFreq)),colorNames(2));
    title(hPlots(5,pos),'H vs M: AI(b) & AO(r)');
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
function y=getAUC(x,condition)

if ~exist('condition','var');       condition = 'A';                    end

if (size(x,1)~=2) || (size(x,2)~=4)
    error('Data matrix has inconsistent size');
end
xSize = numel(size(x{1,1})); % 2 for FR, 3 for LFP
y = cell(2,2);

for i=1:2
    for j=1:2
        if strcmp(condition,'A')
            d1 = abs(x{i,j}); d2 = abs(x{i,j+2});
        elseif strcmp(condition,'P')
            d1 = angle(x{i,j}); d2 = angle(x{i,j+2});
        end
        numElecs = size(d1,1);
        if size(d2,1)~=numElecs
            error('Number of electrodes are different');
        end
        
        if xSize==2 % Firing Rate
            tmpData = zeros(numElecs,1);
            for k=1:numElecs
                tmpData(k) = getDPrime(d1(k,:),d2(k,:));
                %tmpData(k) = ROCAnalysis(d1(k,:),d2(k,:));
            end
        else % LFP power or phase
            numFreqs = size(d1,2);
            tmpData = zeros(numElecs,numFreqs);
            for k=1:numElecs
                for f=1:numFreqs
                    tmpData(k,f) = getDPrime(squeeze(d1(k,f,:)),squeeze(d2(k,f,:)));
                    %tmpData(k,f) = ROCAnalysis(squeeze(d1(k,f,:)),squeeze(d2(k,f,:)));
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