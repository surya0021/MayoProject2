% This program can be used to run saveSegmentedDataMayoVer2 for all the
% sessions

[fileNameStringList,monkeyNameList] = getAttentionExperimentDetails;

folderSourceString='E:\Mayo';    %'C:\Supratim\Projects\MayoProject\';
% instructionTrialFlag=0;

for m=1:length(monkeyNameList)
    
    for i=1:length(fileNameStringList{m})
        saveSegmentedDataMayoVer2(fileNameStringList{m}{i},folderSourceString);
    end
end