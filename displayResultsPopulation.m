% conditionType: 'V', 'N' or 'I' (valid, neutral or invalid)
% targetOnsetMatchingChoice: 1 - nothing, 2 - numtrials, 3 - mean matching (default)
% numTrialCutoff - only select sessions which have more than these number of trials
% TWNum - TW product for Multi-taper analysis

% Modified from displayResultsSingleElectrode. This program shows results
% after first combining data from all electrodes recorded from the same
% session. Functions are taken from displayHvsMPopulation2 used in the
% first project (MayoProject repository; published in Prakash et al., 2021,
% Cerebral Cortex)

% Options for combining data are as follows
% transformType - 1: simple averaging across sides, 2: Mean-difference, 3: LDA-uncorrelated, 4: LDA-covariance
% numFolds: Folds for cross validation. Set numFolds to 1 if you do not want to do cross-validation
% useEqualStimRepsFlag: set to 1 to equalize the stimulus repeats for the two classes
% numElectrodesToUse: Set as empty to use all electrodes. This option allows us to work with a smaller set of electrodes
% regFlag - type of LDA regularization

function displayResultsPopulation(conditionType,targetOnsetMatchingChoice,numTrialCutoff,TWNum,transformType,colorToUse,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag)

if ~exist('targetOnsetMatchingChoice','var'); targetOnsetMatchingChoice=3; end
if ~exist('numTrialCutoff','var');            numTrialCutoff=10;        end
if ~exist('TWNum','var');                   TWNum=3;                    end
if ~exist('transformType','var');           transformType=1;            end
if ~exist('colorToUse','var');              colorToUse='r';             end
if ~exist('numFolds','var');                numFolds=5;                 end
if ~exist('useEqualStimRepsFlag','var');    useEqualStimRepsFlag=0;     end
if ~exist('numElectrodesToUse','var');      numElectrodesToUse=[];      end
if ~exist('regFlag','var');                 regFlag=0;                  end

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Condition Indices %%%%%%%%%%%%%%%%%%%%%%%%%%
% The order of the 12 conditions is as follows: {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'};
if strcmpi(conditionType(1),'V')
    conditionsToUse = [1 2 5 6];
elseif strcmpi(conditionType(1),'N')
    conditionsToUse = 9:12;
elseif strcmpi(conditionType(1),'I')
    conditionsToUse = [3 4 7 8];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The conditions are H0, H1, M0 and M1. Here 0 - attend left and 1 - right.
% We compare the following (same convention as used for single electrode analysis

indicesToCompare{1,1} = [1 2]; % L vs R for Hits
indicesToCompare{1,2} = [3 4]; % L vs R for Misses
indicesToCompare{2,1} = [1 3]; % H vs M for Attend L
indicesToCompare{2,2} = [2 4]; % H vs M for Attend R

% colorsForComparison{1,1} = 'k';
% colorsForComparison{1,2} = [0.5 0.5 0.5]; % Gray
% colorsForComparison{2,1} = 'b';
% colorsForComparison{2,2} = 'r';

legendForComparison{1,1} = 'H(LvR)';
legendForComparison{1,2} = 'M(LvR)';
legendForComparison{2,1} = 'L(HvM)';
legendForComparison{2,2} = 'R(HvM)';

hPlots = getPlotHandles(2,2,[0.05 0.05 0.9 0.9],0.05,0.1,0);

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
        for cType2=1:2
            tmpFiringRates{cType2,k} = allFiringRates0{i}{cType2,k}(:,goodStimNums{i}{k});
            tmpMTPower{cType2,k} = allMTPower0{i}{cType2,k}(:,:,goodStimNums{i}{k});
        end
        tmpAllNumTrials(k) = length(goodStimNums{i}{k});
        tmpTargetOnsetTimes{k} = allTargetOnsetTimes0{i}{k}(goodStimNums{i}{k});
    end
    allFiringRates{i} = tmpFiringRates;
    allMTPower{i} = tmpMTPower;
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

goodFiringRates = allFiringRates(goodSessionList);
goodMTPower = allMTPower(goodSessionList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Transform Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data from all electrodes in a session need to be combined by multiplying
% by a suitable weight vector. The weight vector is chosen to maximise the
% distance between two classes (Attend-Left vs Attend-Out or Hit vs Miss)
freqPosToUse = 1:5:length(freqValsMT);
freqValsMTToUse = freqValsMT(freqPosToUse);
numFreqs = length(freqPosToUse);
dPrimeFRAll = zeros(2,2,numGoodSessions);
dPrimeMTAll = zeros(2,2,numGoodSessions,numFreqs);

for cType1 = 1:2
    for cType2=1:2
        tmpConditions = indicesToCompare{cType1,cType2};
        disp(['Working on ' legendForComparison{cType1,cType2}]);

        for s=1:numGoodSessions    
            disp(['Working on session: ' num2str(s)]);
            
            % Firing Rate
            frData1 = goodFiringRates{s}(:,conditionsToUse(tmpConditions(1)));
            frData2 = goodFiringRates{s}(:,conditionsToUse(tmpConditions(2)));          
            dPrimeFRAll(cType1,cType2,s) = transformData(frData1,frData2,transformType,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag);
            
            % MT Power
            mtData1 = goodMTPower{s}(:,conditionsToUse(tmpConditions(1)));
            mtData2 = goodMTPower{s}(:,conditionsToUse(tmpConditions(2)));
            for k=1:numFreqs
                fPos=freqPosToUse(k);
                mtData1tmp{1} = squeeze(mtData1{1}(:,fPos,:)); mtData1tmp{2} = squeeze(mtData1{2}(:,fPos,:));
                mtData2tmp{1} = squeeze(mtData2{1}(:,fPos,:)); mtData2tmp{2} = squeeze(mtData2{2}(:,fPos,:));
                dPrimeMTAll(cType1,cType2,s,k) = transformData(mtData1tmp,mtData2tmp,transformType,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag); 
            end
        end
     
        plotData(hPlots(cType1,cType2),freqValsMTToUse,squeeze(dPrimeMTAll(cType1,cType2,:,:)),colorToUse); % MT Power
        plotData(hPlots(cType1,cType2),freqValsMTToUse,repmat(squeeze(dPrimeFRAll(cType1,cType2,:)),1,numFreqs),colorToUse); % Firing Rate
        title(hPlots(cType1,cType2),legendForComparison{cType1,cType2});
    end
end
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

% Functions for Population Analysis taken from previous project
function [dPrime,mD] = transformData(inpData1,inpData2,transformType,numFolds,useEqualStimRepsFlag,numElectrodesToUse,regFlag)

n1=size(inpData1{1},1); n2=size(inpData1{2},1); % Number of electrodes
data1 = [inpData1{1};inpData1{2}]; % Concatenate data from both sides
data2 = [inpData2{1};inpData2{2}];

eList = 1:(n1+n2);
eSideList = [zeros(1,n1) ones(1,n2)];
if ~isempty(numElectrodesToUse)
    if n1+n2>numElectrodesToUse % Use a subset of electrodes
        eIndices = randperm(n1+n2);
        eIndices = sort(eIndices(1:numElectrodesToUse));
        eList = eList(eIndices);
        n1 = length(find(eSideList(eIndices)==0));
        n2 = length(find(eSideList(eIndices)==1));
    end
end

data1 = data1(eList,:);
data2 = data2(eList,:);

[testingIndices1,testingIndices2] = getIndices(size(data1,2),size(data2,2),numFolds,useEqualStimRepsFlag);

% Get the weightVector and projections of data1 and data2
[p1,p2,~,allIndices1,allIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2,regFlag);
p1 = p1(allIndices1);
p2 = p2(allIndices2);
dPrime = getDPrime(p1,p2);
mD = mean(p1) - mean(p2);
end
function [testingIndices1,testingIndices2] = getIndices(N1,N2,numFolds,useEqualStimRepsFlag)

allIndices1 = randperm(N1);
allIndices2 = randperm(N2);

if useEqualStimRepsFlag
    N1=min(N1,N2);
    N2=N1;
end

allIndices1 = sort(allIndices1(1:N1));
allIndices2 = sort(allIndices2(1:N2));
allIndices = [allIndices1 allIndices2];

testingIndices1 = cell(1,numFolds);
testingIndices2 = cell(1,numFolds);

if numFolds==1 % No cross Validation
    testingIndices1{1} = allIndices1;
    testingIndices2{1} = allIndices2;
else
    Y = [zeros(N1,1) ; ones(N2,1)];
    cvp = cvpartition(Y,'KFold',numFolds); % Partition data
    
    for i=1:numFolds
        testingIDs = find(cvp.test(i)==1);
        testingIndices1{i} = allIndices(testingIDs(testingIDs<=N1));
        testingIndices2{i} = allIndices(testingIDs(testingIDs>N1));
    end
end
end
function [projections1,projections2,weightVector,fullSetIndices1,fullSetIndices2] = getProjectionsAndWeightVector(data1,data2,testingIndices1,testingIndices2,transformType,n1,n2,regFlag)

numFolds = length(testingIndices1);
fullSetIndices1=[]; fullSetIndices2=[];
for i=1:numFolds
    fullSetIndices1 = cat(2,fullSetIndices1,testingIndices1{i});
    fullSetIndices2 = cat(2,fullSetIndices2,testingIndices2{i});
end
fullSetIndices1 = sort(fullSetIndices1);
fullSetIndices2 = sort(fullSetIndices2);

projections1 = zeros(size(data1,2),1);
projections2 = zeros(size(data2,2),1);

D = size(data1,1);
weightVectorTMP = zeros(D,numFolds);

for i=1:numFolds
    t1 = testingIndices1{i};
    t2 = testingIndices2{i};
    
    if transformType==1 % Simple Averaging of sides
        weightVectorTMP(:,i) = [ones(n1,1)/n1; -ones(n2,1)/n2];
        
    else
        if numFolds==1 % No cross validation. Train and test on the same data
            train1 = t1;
            train2 = t2;
        else
            train1 = setdiff(fullSetIndices1,t1);
            train2 = setdiff(fullSetIndices2,t2);
        end
        d1 = data1(:,train1); d2 = data2(:,train2);
        m1 = size(d1,2); m2 = size(d2,2);
        
        if transformType==2 % Only mean difference
            weightVectorTMP(:,i) = mean(d1,2) - mean(d2,2);
            
        elseif transformType==3 % LDA for uncorrelated case       
            meanDiff = mean(d1,2) - mean(d2,2);
            var1 = var(d1,[],2);
            var2 = var(d2,[],2);
            pooledVar = ((m1-1)*var1 + (m2-1)*var2)/(m1+m2-2); %Pooled Variance
            
            weightVectorTMP(:,i) = meanDiff./pooledVar;
            
        elseif transformType==4 % LDA

            label1 = repmat({'Class1'},m1,1);
            label2 = repmat({'Class2'},m2,1);
            labelAll = cat(1,label1,label2);
            dataAll = cat(2,d1,d2);

            if regFlag==0 % No regularization
                Mdl = fitcdiscr(dataAll',labelAll);
                
            elseif regFlag==1 % Only optimize gamma
                Mdl = fitcdiscr(dataAll',labelAll);
                [err,gamma,~,~] = cvshrink(Mdl,'NumGamma',20);
                [~,minIndex] = min(err);
                if minIndex>1
                    Mdl.Gamma = gamma(minIndex); % Changing gamma changes the model weights also
                end
                
            elseif regFlag==2 % Optimize gamma and delta
                Mdl = fitcdiscr(dataAll',labelAll);
                [err,gamma,delta,~] = cvshrink(Mdl,'NumGamma',20,'numDelta',20);
                minerr = min(err(:));
                [x,y] = find(err==minerr);
                if x(1)>1 || y(1)>1
                    Mdl.Gamma = gamma(x(1)); % Take the smallest gamma
                    Mdl.Delta = delta(x(1),y(1)); % Take the smallest delta
                end
                
            elseif regFlag==3 % Hyper optimize gamma and delta
                rng(1)
                myOpts.AcquisitionFunctionName = 'expected-improvement-plus';
                myOpts.ShowPlots = 0;
                myOpts.Verbose = 0;
                Mdl = fitcdiscr(X,Y,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',myOpts);
            end
            
            weightVectorTMP(:,i) = Mdl.Coeffs(1,2).Linear;
        end
    end
    projections1(t1) = data1(:,t1)' *  weightVectorTMP(:,i);
    projections2(t2) = data2(:,t2)' *  weightVectorTMP(:,i);
end

weightVector = mean(weightVectorTMP,2);
end