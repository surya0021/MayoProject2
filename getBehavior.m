function [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg,reactionTimeMS] = getBehavior(fileNameString,folderSourceString,neutralTrialFlag)

if ~exist('folderSourceString','var');   folderSourceString='G:\Mayo';       end
if ~exist('neutralTrialFlag','var');    neutralTrialFlag = 1;   end
folderNameIn = fullfile(folderSourceString,'Data','extractedData');

fileNameCDS = [fileNameString '_Codes'];
CDS = load(fullfile(folderNameIn,fileNameCDS));
CDS = CDS.(fileNameCDS);

fileNameDAT = strcat(fileNameString, '_DAT');
DAT = load(fullfile(folderNameIn,fileNameDAT));
DAT = DAT.(fileNameDAT);

goodIndexList = getGoodIndices(CDS,DAT,[],neutralTrialFlag);
[~,~,~,orientationChangeDeg,~,reactionTimeMS] = getInfoDATFile(DAT);

uniqueOrientationChangeDeg = unique(orientationChangeDeg);

numConditions = length(goodIndexList);
numOrientations = length(uniqueOrientationChangeDeg);

responseMatrix = zeros(numConditions,numOrientations);
for i=1:numConditions
    x = orientationChangeDeg(goodIndexList{i});
    for j=1:numOrientations
        responseMatrix(i,j) = length(find(x==uniqueOrientationChangeDeg(j)));
    end
end

if neutralTrialFlag
    perCorrect = zeros(6,numOrientations); % order: 0V, 1V, 0I, 1I, 0N, 1N
    for i=1:4
        perCorrect(i,:) = responseMatrix(i,:) ./ (responseMatrix(i,:)+responseMatrix(i+4,:));
    end
    perCorrect(5,:) = responseMatrix(9,:) ./ (responseMatrix(9,:)+responseMatrix(11,:));
    perCorrect(6,:) = responseMatrix(10,:) ./ (responseMatrix(10,:)+responseMatrix(12,:));
else
    perCorrect = zeros(5,numOrientations); % order: 0V, 1V, 0I, 1I, N
    for i=1:4
        perCorrect(i,:) = responseMatrix(i,:) ./ (responseMatrix(i,:)+responseMatrix(i+4,:));
    end
    perCorrect(5,:) = responseMatrix(9,:) ./ (responseMatrix(9,:)+responseMatrix(10,:));
end
end