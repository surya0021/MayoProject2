% Created December 5, 2012.  JPM
% extracts and groups trials of
% various conditions of interest, including instruction trials, correct,
% etc.
% indKeptTrials is logical index of CODES input

% From getSpikesCodes: CDS/codes format = TrialStart FixOn FixateOn CueOn StimOn StimOff Sacc TrialEnd ;

function [indKeptTrials, indCuedLoc0, indCuedLoc1,indUncued] = getTrialTypes (CDS, correct, miss, instruct)
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
