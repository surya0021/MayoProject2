% This program is taken from MayoProjectPrograms. Several other functions
% that this program calls were separate files in the previous project, but
% are just concatenated in this program.

% Option added to save this data locally. That way, the same program can be
% used even when the main Data folder is unavailable.

function [perCorrect,uniqueOrientationChangeDeg,goodIndexList,orientationChangeDeg,reactionTimeMS,targetOnTimeMS] = getBehavior(fileNameString,folderSourceString,neutralTrialFlag)

if ~exist('folderSourceString','var');   folderSourceString='G:\Mayo';  end
if ~exist('neutralTrialFlag','var');    neutralTrialFlag = 1;           end

%%%%%%%%%%%%% check whether file exists within a local folder %%%%%%%%%%%%%
savedDataFolder = fullfile(pwd,'savedData','behaviorData');
savedDataFile = fullfile(savedDataFolder,[fileNameString '.mat']);

if exist(savedDataFile,'file')
    load(savedDataFile,'perCorrect','uniqueOrientationChangeDeg','goodIndexList','orientationChangeDeg','reactionTimeMS','targetOnTimeMS');
else
    
    folderNameIn = fullfile(folderSourceString,'Data','extractedData');
    
    fileNameCDS = [fileNameString '_Codes'];
    CDS = load(fullfile(folderNameIn,fileNameCDS));
    CDS = CDS.(fileNameCDS);
    
    fileNameDAT = strcat(fileNameString, '_DAT');
    DAT = load(fullfile(folderNameIn,fileNameDAT));
    DAT = DAT.(fileNameDAT);
    
    goodIndexList = getGoodIndices(CDS,DAT,[],neutralTrialFlag);
    [~,~,targetOnTimeMS,orientationChangeDeg,~,reactionTimeMS] = getInfoDATFile(DAT);
    
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
    
    % Save data
    makeDirectory(savedDataFolder);
    save(savedDataFile,'perCorrect','uniqueOrientationChangeDeg','goodIndexList','orientationChangeDeg','reactionTimeMS','targetOnTimeMS');
end
end
function goodIndexList = getGoodIndices(CDS,DAT,instructionTrialFlag,neutralTrialFlag)

% This program is directly copied from MayoProjectPrograms 
% instructionTrialFlag - 1 generates the indices for intructionTrials
% neutralTrialFlag - 1 splits the neutral condition and generates indices for left and right change(target)
% neutral trials

if ~exist('instructionTrialFlag','var') || isempty(instructionTrialFlag) ;    instructionTrialFlag=0;     end
if ~exist('neutralTrialFlag','var');    neutralTrialFlag=0;     end

[DAT2,CueOnCode] = getInfoDATFile(DAT);
DAT2 = DAT2(1:length(CDS)); % DAT file sometimes has more entries than the CDS file. The last block for attLoc1 has not been used for analysis.
CueOnCode = CueOnCode(1:length(CDS));

%trialEndCode = cell2mat(cellfun(@(x) x.trialEnd.data, DAT2,'UniformOutput',false));
trialEndCode = cell2mat(cellfun(@(x) x{8},CDS(:,2),'UniformOutput',false))';

isValidTrial = cell2mat(cellfun(@(x) x.trial.data.validTrial,DAT2,'UniformOutput',false));
% isInstructTrial = cell2mat(cellfun(@(x) x.trial.data.instructTrial,DAT2,'UniformOutput',false));
isInstructTrial = (CueOnCode>0); % Sometimes a cue is generated even if the trial is not an instruction trial. Therefore, any trial which contains the cueOn field is considered an instruction trial
attendLoc = cell2mat(cellfun(@(x) x.trial.data.attendLoc,DAT2,'UniformOutput',false));
correctLoc=cell2mat(cellfun(@(x) x.trial.data.correctLoc,DAT2,'UniformOutput',false));
isCatchTrial = cell2mat(cellfun(@(x) x.trial.data.catchTrial,DAT2,'UniformOutput',false));

% Hit Trials
hitIndices = (trialEndCode==0) & (isInstructTrial==instructionTrialFlag) & (isCatchTrial==0);
goodIndexList{1} = find(hitIndices & (isValidTrial==1) & (attendLoc==0)); % H0V
goodIndexList{2} = find(hitIndices & (isValidTrial==1) & (attendLoc==1)); % H1V
goodIndexList{3} = find(hitIndices & (isValidTrial==0) & (attendLoc==0)); % H0I
goodIndexList{4} = find(hitIndices & (isValidTrial==0) & (attendLoc==1)); % H1I

% Miss Trials
missIndices = (trialEndCode==2) & (isInstructTrial==instructionTrialFlag) & (isCatchTrial==0);
goodIndexList{5} = find(missIndices & (isValidTrial==1) & (attendLoc==0)); % M0V
goodIndexList{6} = find(missIndices & (isValidTrial==1) & (attendLoc==1)); % M1V
goodIndexList{7} = find(missIndices & (isValidTrial==0) & (attendLoc==0)); % M0I
goodIndexList{8} = find(missIndices & (isValidTrial==0) & (attendLoc==1)); % M1I

if neutralTrialFlag
    goodIndexList{9} = find(hitIndices & (isValidTrial==1) & (attendLoc==2) & (correctLoc==0)); % H0N
    goodIndexList{10} = find(hitIndices & (isValidTrial==1) & (attendLoc==2) & (correctLoc==1)); % H1N
    goodIndexList{11} = find(missIndices & (isValidTrial==1) & (attendLoc==2) & (correctLoc==0)); % M0N
    goodIndexList{12} = find(missIndices & (isValidTrial==1) & (attendLoc==2) & (correctLoc==1)); % M1N
else
    goodIndexList{9} = find(hitIndices & (isValidTrial==1) & (attendLoc==2)); % HN
    goodIndexList{10} = find(missIndices & (isValidTrial==1) & (attendLoc==2)); % MN
end

% Doing the same thing using Patrick's code
[~, indHitLoc0, indHitLoc1,indHitNeutral] = getTrialTypes (CDS,1,0,0);
catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitLoc0 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indHitLoc0 ); % logical index of valid trials for Hit Loc 0 trials
goodIndexList2{1} = (indHitLoc0(validList & ~catchList))'; % HOV
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indHitLoc0 ); % logical index of invalid trials for Hit Loc 0 trials
goodIndexList2{3} = (indHitLoc0(invalidList & ~catchList))'; % HOI

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitLoc1 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indHitLoc1 ); % logical index of valid trials for Hit Loc 1 trials
goodIndexList2{2} = (indHitLoc1(validList & ~catchList))'; % H1V
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indHitLoc1 ); % logical index of invalid trials for Hit Loc 1 trials
goodIndexList2{4} = (indHitLoc1(invalidList & ~catchList))'; % H1I

[~, indMissLoc0, indMissLoc1,indMissNeutral] = getTrialTypes (CDS,0,1,0);
catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissLoc0 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indMissLoc0 ); % logical index of valid trials for Miss Loc 0 trials
goodIndexList2{5} = (indMissLoc0(validList & ~catchList))'; % MOV
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indMissLoc0 ); % logical index of invalid trials for Miss Loc 0 trials
goodIndexList2{7} = (indMissLoc0(invalidList & ~catchList))'; % MOI

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissLoc1 ); % logical index of catch trials with no stimulus change
validList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==1, indMissLoc1 ); % logical index of valid trials for Miss Loc 1 trials
goodIndexList2{6} = (indMissLoc1(validList & ~catchList))'; % M1V
invalidList= arrayfun(@(x) DAT2{x}.trial.data.validTrial(1)==0, indMissLoc1 ); % logical index of invalid trials for Miss Loc 1 trials
goodIndexList2{8} = (indMissLoc1(invalidList & ~catchList))'; % M1I

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indHitNeutral ); % logical index of catch trials with no stimulus change
goodIndexList2{9} = (indHitNeutral(~catchList))'; % HN

catchList= arrayfun(@(x) DAT2{x}.trial.data.catchTrial(1)==1, indMissNeutral ); % logical index of catch trials with no stimulus change
goodIndexList2{10} = (indMissNeutral(~catchList))'; % MN

if ~instructionTrialFlag && ~neutralTrialFlag
    if ~isequal(goodIndexList,goodIndexList2)
        error('Index Lists do not match...');
    end
end
end
function [DAT2,CueOnCode,targetOnTimeMS,orientationChangeDeg,saccadeTimeMS,reactionTimeMS] = getInfoDATFile(DAT)
i=1; continueFlag=1;
while(continueFlag)
    tnum = (['t', num2str(i)]);
    if isfield(DAT,tnum)
        DAT2{i} = DAT.(tnum); %#ok<*AGROW>
        if isfield(DAT2{i},'cueOn')
            CueOnCode(i) = 1;
        else
            CueOnCode(i) = 0;
        end
        
        if isfield(DAT2{i}.trial.data,'targetOnTimeMS')
            targetOnTimeMS(i) = DAT2{i}.trial.data.targetOnTimeMS; % From stimulus onset
            orientationChangeDeg(i) = DAT2{i}.trial.data.change;
        else
           targetOnTimeMS(i) = 0;
           orientationChangeDeg(i) = 0;
        end
        
        if isfield(DAT2{i},'saccade')
            saccadeTimeMS(i) = DAT2{i}.saccade.timeMS - DAT2{i}.trialStart.timeMS;
            if isfield(DAT2{i},'target')
                reactionTimeMS(i) = DAT2{i}.saccade.timeMS - DAT2{i}.target.timeMS(1);
            else
                reactionTimeMS(i) = 0;
            end
        else
           saccadeTimeMS(i) = 0;
           reactionTimeMS(i) = 0;
        end
        
        i=i+1;
    else
        continueFlag=0;
    end
end
end
function [indKeptTrials, indCuedLoc0, indCuedLoc1,indUncued] = getTrialTypes (CDS, correct, miss, instruct)

% Created December 5, 2012.  JPM
% extracts and groups trials of
% various conditions of interest, including instruction trials, correct,
% etc.
% indKeptTrials is logical index of CODES input

% From getSpikesCodes: CDS/codes format = TrialStart FixOn FixateOn CueOn StimOn StimOff Sacc TrialEnd ;

% NOTE indKeptTrials is a LOGICAL INDEX, other indices are position indices

nTrials = size(CDS,1); % number of total trials from data file

% Extract the labels from the Codes
CueOnCode = cell2mat( cellfun( @(x) x{4}, CDS(:,2), 'UniformOutput', false ) ); % Not a logical index. Extracts values from original CODES (eg, 0, 1, 2, NaN) for all trials
TrialEndCode = cell2mat( cellfun( @(x) x{8}, CDS(:,2), 'UniformOutput', false ) ); % trial end codes (0: correct; 1; 2: miss ,3 or 4)


% extract only certain types of trials (the ones with "true" value)

if correct
    indCorrect = TrialEndCode == 0; % gives LOGICAL index of all trials, 1 when true and 0 else
else
    indCorrect = zeros(nTrials,1);
end

if miss
    indMiss = TrialEndCode == 2; % gives LOGICAL index of all trials, 1 when true and 0 else
else
    indMiss = zeros(nTrials,1);
end

if instruct
    indInstr = ~isnan(CueOnCode); % gives LOGICAL index of all trials, al1 when its an instruction tris (is not NaN) and 0 when NaN, not an instrut trials
else
    indInstr = isnan(CueOnCode);
end

indKeptTrials = (indCorrect | indMiss) & indInstr ; % logical indexes of the kept trials


b = find( ~isnan(CueOnCode)); % fill in first cued location value with cue in first trial
CuedLoc = CueOnCode(b(1));
cuedloc = CuedLoc;

% CueOnCode = 0, 1, 2 or NaN
for g=2:length(CueOnCode) % go thru all cue codes, starting with trial 2, and fill in Cued Location
     if ~isnan(CueOnCode(g)) && CueOnCode(g) ~= cuedloc
         cuedloc = CueOnCode(g);
     end
     CuedLoc = [CuedLoc; cuedloc];
end

% we need the cued locations of the kept trials only

indCuedLoc0 = find (indKeptTrials & CuedLoc == 0); % All kept trials (instruction and normal) where stimulus is presented at Location 0, usually left
indCuedLoc1 = find (indKeptTrials & CuedLoc == 1); % All kept trials (instruction and normal) where stimulus is presented at Location 1, usually right
indUncued  =  find (indKeptTrials & CuedLoc == 2); % Uncued

end