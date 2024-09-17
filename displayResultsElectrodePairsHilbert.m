% This program computes and displays the mean difference and d-prime of shuffle corrected LFP power correlation
% and PPC calculated using Hilbert transform between attention and behavioral conditions.

function [allMeansPPCHT,allMeansCorrHTCorrectedAndRaw,allDPrimesPPCHT,allDPrimesCorrHTCorrectedAndRaw,freqValsHT] = displayResultsElectrodePairsHilbert(conditionType,freqHalfBandwidth,targetOnsetMatchingChoice,numTrialCutoff,displayFlag)
if ~exist('freqHalfBandwidth','var');            freqHalfBandwidth = 5;       end
if ~exist('targetOnsetMatchingChoice','var');    targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');               numTrialCutoff=10;           end
if ~exist('displayFlag','var');                  displayFlag = 1;             end


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


%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[allPPCHT0,allPowerCorrHT0,~,freqValsHT,~,~,allTargetOnsetTimes0] = getHilbertPairwiseMeasures(freqHalfBandwidth);

numSessions = length(allTargetOnsetTimes0);
numConditions = length(allTargetOnsetTimes0{1});

%%%%%%%%%%%%%%%%%%%% Get good StimIndices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetTimeBinWidthMS = 250;
goodStimNums = getGoodStimNums(allTargetOnsetTimes0,targetOnsetMatchingChoice,targetTimeBinWidthMS);

%%%%%%%%%%%%%%%%%%%% Select only good stimIndices %%%%%%%%%%%%%%%%%%%%%%%%%
allPPCHT = cell(2,numConditions,numSessions);
allPowerCorrHT = cell(2,numConditions,numSessions);
allNumTrials = cell(1,numSessions);
allTargetOnsetTimes = cell(1,numSessions);

for i=1:numSessions
    tmpAllNumTrials = zeros(1,numConditions);
    tmpTargetOnsetTimes = cell(1,numConditions);

    for k=1:numConditions
        for j=1:2
            allPPCHT{j,k,i} = allPPCHT0{j,k,i}(goodStimNums{i}{k},:,:);
            allPowerCorrHT{j,k,i} = allPowerCorrHT0{j,k,i}(goodStimNums{i}{k},:,:,:);
        end
        tmpAllNumTrials(k) = length(goodStimNums{i}{k});
        tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
    end
    allNumTrials{i} = tmpAllNumTrials;
    allTargetOnsetTimes{i} = tmpTargetOnsetTimes;
end


%%%%%%%%%%%%%%%%%%%% Select sessions with numTrials>=cutoff %%%%%%%%%%%%%%%
allNumTrialsMatrix = cell2mat(allNumTrials');
minTrialsConditions = min(allNumTrialsMatrix(:,conditionsToUse),[],2);

badSessionList = find(minTrialsConditions<=numTrialCutoff);
goodSessionList = setdiff(1:length(allNumTrials),badSessionList);
numGoodSessions = length(goodSessionList);
disp(['Discarded sessions: ' num2str(badSessionList')]);

goodPPCHT = allPPCHT(:,:,goodSessionList);
goodCorrHT = allPowerCorrHT(:,:,goodSessionList);

%%%%%%%%%%%%%%%%%%%% Get Means and DPrimes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pairwise data is for conditions H0, H1, M0 and M1. This needs to be
% converted to Hit(InVsOut), Miss(InVsOut), In(HitVsMiss) and Out(HitVsMiss)
allMeansPPCHT = [];    allDPrimesPPCHT = [];
allMeansCorrHT = [];    allDPrimesCorrHT = [];

nType=2; % 1:PPCHT  2:CorrHT

for dataType=1:nType
    clear allMeans allDPrimes pairwiseMeasureData
    if dataType==1
        pairwiseMeasureData = goodPPCHT(:,conditionsToUse,:);
        allMeans0 = cell(2,4);
        allDPrimes0 = cell(2,4);
        for side=1:2
            for c=1:4 % Comparison
                tmpIndexToCompare = indicesToCompare{c};
                disp(['side' num2str(side) '_' originalConditionsList{tmpIndexToCompare(1)} '-' originalConditionsList{tmpIndexToCompare(2)}]);

                meanDataAllSessions=[];
                dPrimeDataAllSessions=[];
                for i=1:numGoodSessions
                    tmpData1 = pairwiseMeasureData{side,tmpIndexToCompare(1),i};
                    tmpData2 = pairwiseMeasureData{side,tmpIndexToCompare(2),i};

                    [~,numFreqs,numPairs] = size(tmpData1);

                    clear tmpMeanData tmpDPrime
                    tmpMeanData = zeros(numFreqs,numPairs);
                    tmpDPrimeData = zeros(numFreqs,numPairs);
                    for f=1:numFreqs
                        for p=1:numPairs
                            x1 = squeeze(tmpData1(:,f,p));
                            x2 = squeeze(tmpData2(:,f,p));
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
        allMeansPPCHT = cat(3,allMeans{:});
        allDPrimesPPCHT = cat(3,allDPrimes{:});
    elseif dataType==2
        pairwiseMeasureData = goodCorrHT(:,conditionsToUse,:);
        allMeans0 = cell(2,4);
        allDPrimes0 = cell(2,4);
        for side=1:2
            for c=1:4 % Comparison
                tmpIndexToCompare = indicesToCompare{c};
                disp(['side' num2str(side) '_' originalConditionsList{tmpIndexToCompare(1)} '-' originalConditionsList{tmpIndexToCompare(2)}]);

                meanDataAllSessions=[];
                dPrimeDataAllSessions=[];
                for i=1:numGoodSessions
                    tmpData1 = pairwiseMeasureData{side,tmpIndexToCompare(1),i};
                    tmpData2 = pairwiseMeasureData{side,tmpIndexToCompare(2),i};

                    [~,numFreqs,numPairs,numCorrTypes] = size(tmpData1);

                    clear tmpMeanData tmpDPrime
                    tmpMeanData = zeros(numFreqs,numPairs,numCorrTypes);
                    tmpDPrimeData = zeros(numFreqs,numPairs,numCorrTypes);
                    for f=1:numFreqs
                        for p=1:numPairs
                            x1 = squeeze(tmpData1(:,f,p,:));
                            x2 = squeeze(tmpData2(:,f,p,:));
                            tmpMeanData(f,p,:) = mean(x1,'omitnan')-mean(x2,'omitnan');
                            tmpDPrimeData(f,p,:) = getDPrime(x1,x2);
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
        allMeansCorrHT = cat(4,allMeans{:});
        allDPrimesCorrHT = cat(4,allDPrimes{:});
        
        allMeansCorrHTCorrectedAndRaw = allMeansCorrHT(:,:,1:2,:);
        allDPrimesCorrHTCorrectedAndRaw = allDPrimesCorrHT(:,:,1:2,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayFlag
    figure
    hPlots = getPlotHandles(2,4,[0.05 0.05 0.9 0.9],0.05,0.1,0);
    %%%%%%%%%%%%%%%%%%%%%%%% Plotting correlation %%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        %plot mean change in firing rate correlation with power only
        pCorrHT(i,1) = plotData(hPlots(1,2*plotPos-1),freqValsHT,allMeansCorrHT(:,:,1,i)',colorsForComparison{i});
        pCorrHT(i,2) = plotData(hPlots(1,2*plotPos),freqValsHT,allDPrimesCorrHT(:,:,1,i)',colorsForComparison{i});
    end
    
    for i=1:2
        legend(hPlots(1,i),pCorrHT(1:2,2),legendForComparison(1:2),'autoupdate','off');
        legend(hPlots(1,i+2),pCorrHT(3:4,2),legendForComparison(3:4),'autoupdate','off');
        plot(hPlots(1,(2*i-1)),freqValsHT,zeros(1,length(freqValsHT)),'k--','linewidth',1.5);
        plot(hPlots(1,2*i),freqValsHT,zeros(1,length(freqValsHT)),'k--','linewidth',1.5);
        ylabel(hPlots(1,(2*i-1)),'\Delta Correlation');
        ylabel(hPlots(1,2*i),'dPrime');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the PPC %%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:4
        if (i<3); plotPos = 1; else; plotPos=2; end
        pPPC(i,1) = plotData(hPlots(2,2*plotPos-1),freqValsHT,allMeansPPCHT(:,:,i)',colorsForComparison{i},0); %#ok<*AGROW>
        pPPC(i,2) = plotData(hPlots(2,2*plotPos),freqValsHT,allDPrimesPPCHT(:,:,i)',colorsForComparison{i},0);
    end

    for i=1:2
        legend(hPlots(2,i),pPPC(1:2,i),legendForComparison(1:2),'AutoUpdate','off');
        legend(hPlots(2,i+2),pPPC(3:4,i),legendForComparison(3:4),'AutoUpdate','off');
    end
    for i=1:2
        plot(hPlots(2,(2*i-1)),freqValsHT,zeros(1,length(freqValsHT)),'k--');
        plot(hPlots(2,2*i),freqValsHT,zeros(1,length(freqValsHT)),'k--');
        ylabel(hPlots(2,(2*i-1)),'\Delta PPC');
        ylabel(hPlots(2,2*i),'dPrime');
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