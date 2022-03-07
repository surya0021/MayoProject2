% runDisplayBehaviorAndSaveData
% This is the main program that displays the behavioral data and optionally
% also saves the intermediate data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDataFlag=0; % Set to 1 if you also want to save data. Otherwise this program just displays behavioral data

selectedOrientations = [2 3]; % Choose only between these two, since we have invalid conditions for only these condition, and since these were overrepresented in the dataset. Previous papers have done the same
colorList = 'rmgkbc';
legendList = [{'0V'} {'1V'} {'0I'} {'1I'} {'0N'} {'1N'}];
numConditions = length(legendList);

%%%%%%%%%%%%%%%%%%%%%%%%%% Get Experimental Details %%%%%%%%%%%%%%%%%%%%%%%
folderSourceString = fileparts(pwd); % Parent folder
fileNameStringList = getAttentionExperimentDetails;
fileNameStringListAll = cat(2,fileNameStringList{1},fileNameStringList{2});
numSessions = length(fileNameStringListAll);

%%%%%%%%%%%%%%%%%%%%%%%% Get electrodes to save %%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveDataFlag
    x=load(fullfile(pwd,'savedData','electrodeArrayListStimulated.mat')); % saved locally from Data/savedDataSummary
    electrodeArrayList = x.electrodeArrayList;
end

allGoodPos = zeros(numSessions,numConditions);
for i=1:numSessions
    fileNameString = fileNameStringListAll{i};
    disp([num2str(i) ': ' fileNameString]);
    [perCorrect,uniqueOrientationChangeDeg] = getBehavior(fileNameString,folderSourceString,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = subplot(5,5,i);
    for j=1:numConditions
        plot(uniqueOrientationChangeDeg,perCorrect(j,:),'color', colorList(j), 'marker','.'); hold on;
    end
    for j=1:numConditions
        plot(uniqueOrientationChangeDeg,perCorrect(j,:),'color', colorList(j));
    end
    plot(uniqueOrientationChangeDeg,0.5+zeros(1,length(uniqueOrientationChangeDeg)),'k');
    
    title(fileNameStringListAll{i});
    ylim([0 1]);
    xlim([0 uniqueOrientationChangeDeg(end)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Choose Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the index out of 2 and 3 which is closest to 0.5
    usefulPerCorrectDeviation = abs(perCorrect(:,selectedOrientations)-0.5);
    for j=1:numConditions
        x = usefulPerCorrectDeviation(j,:);
        tmp = find(x==min(x));
        allGoodPos(i,j) = selectedOrientations(tmp(1));
        plot(uniqueOrientationChangeDeg(allGoodPos(i,j)),perCorrect(j,allGoodPos(i,j)),'color', colorList(j),'marker','o','markerfacecolor',colorList(j));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveDataFlag
        saveDataForAnalysis(fileNameString,folderSourceString,allGoodPos(i,:),electrodeArrayList{i}); %#ok<*UNRCH>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  Display Summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
legend(legendList,'Location','SouthEast');

cueTypeList = [{'Valid'} {'InValid'} {'Neutral'}];
badPosAll=[];

disp('Sessions with different orientation change in left and right');
for i=1:length(cueTypeList)
    d = allGoodPos(:,2*i) - allGoodPos(:,2*(i-1)+1); 
    badPos = (find(abs(d)>0));
    badPosAll = union(badPosAll,badPos);
    disp([cueTypeList{i} ': ' num2str(badPos') ' (N=' num2str(length(badPos)) ')']);
end
disp(['all sessions: ' num2str(badPosAll') ' (N=' num2str(length(badPosAll)) ')']);
