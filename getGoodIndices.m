% This program is directly copied from MayoProjectPrograms which is
% required by saveSegmentedDataMayoVer2 and
% saveSegmentedNeutralDataMayoVer2

% instructionTrialFlag - 1 generates the indices for intructionTrials
% neutralTrialFlag - 1 splits the neutral condition and generates indices for left and right change(target)
% neutral trials

function goodIndexList = getGoodIndices(CDS,DAT,instructionTrialFlag,neutralTrialFlag)

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