% This program is used to find a subset of stimuli for matching hit and
% miss conditions such that the target onset times are matched in
% histograms of size targetTimeBinWidthMS. See code for details.

% allTargetOnsetTimes - cell array of size numSession x 12, where the 12
% conditions are {'H0V','H1V','H0I','H1I','M0V','M1V','M0I','M1I','H0N','H1N','M0N','M1N'}
% targetOnsetMatchingChoice: 1 - no matching, 2 - match num trials, 3 - match mean target onset times

function allGoodStimNums = getGoodStimNums(allTargetOnsetTimes,targetOnsetMatchingChoice,targetTimeBinWidthMS)

numSessions = length(allTargetOnsetTimes);

% Create a 2D matrix with matching hit and miss conditions
allConditionNumbers = [1 2 3 4 9 10; 5 6 7 8 11 12];
numConditions = size(allConditionNumbers,2);

targetOnsetEdges = 500:targetTimeBinWidthMS:5500; % min to max target onset times
numTargetOnsetBins = length(targetOnsetEdges)-1;

allGoodStimNums = cell(1,numSessions);
for i=1:numSessions
    
    goodStimNums = cell(1,2*numConditions);
    for j=1:numConditions
        conditionPos = allConditionNumbers(:,j);
        targetOnsetTimesToUse = allTargetOnsetTimes{i}(conditionPos); % Hit and Miss conditions
        
        if targetOnsetMatchingChoice==1 % Take all trial numbers
            goodStimNums{conditionPos(1)} = 1:length(targetOnsetTimesToUse{1});
            goodStimNums{conditionPos(2)} = 1:length(targetOnsetTimesToUse{2});
            
        elseif targetOnsetMatchingChoice==2 % Take equal number of stimulus repeats for both conditions
            stimNums1 = 1:length(targetOnsetTimesToUse{1});
            stimNums2 = 1:length(targetOnsetTimesToUse{2});
            [goodStimNums{conditionPos(1)},goodStimNums{conditionPos(2)}] = getEqualNumOfIndices(stimNums1,stimNums2);
            
        elseif targetOnsetMatchingChoice==3 % Match the histograms of target onset times
            
            tmpPos1 = []; tmpPos2 = [];
            for k=1:numTargetOnsetBins
                stimNums1 = intersect(find(targetOnsetTimesToUse{1}>=targetOnsetEdges(k)),find(targetOnsetTimesToUse{1}<targetOnsetEdges(k+1)));
                stimNums2 = intersect(find(targetOnsetTimesToUse{2}>=targetOnsetEdges(k)),find(targetOnsetTimesToUse{2}<targetOnsetEdges(k+1)));
                [tmp1,tmp2] = getEqualNumOfIndices(stimNums1,stimNums2);
                tmpPos1 = cat(2,tmp1,tmpPos1);
                tmpPos2 = cat(2,tmp2,tmpPos2);
            end
            goodStimNums{conditionPos(1)} = sort(tmpPos1);
            goodStimNums{conditionPos(2)} = sort(tmpPos2);
        end
    end
    allGoodStimNums{i}=goodStimNums;
end
end
function [x2,y2] = getEqualNumOfIndices(x1,y1)

N1 = length(x1);
N2 = length(y1);

if (N1==0) || (N2==0) % one of the two is an empty array
    x2=[]; y2=[];
    
elseif N1==N2
    x2=x1; y2=y1;
    
elseif N1<N2
    x2 = x1;
    randVals = randperm(N2);
    y2 = y1(sort(randVals(1:N1)));
    
else %N1>N2
    y2 = y1;
    randVals = randperm(N1);
    x2 = x1(sort(randVals(1:N2)));
end
end