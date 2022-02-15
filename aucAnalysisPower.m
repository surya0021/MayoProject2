% This program gets the LFP power and firing rate and ROC analysis data of LFP power and firing rate for single electrode.
% This program is called in plotAucAnalysisPower.m

%parameters to set : 
function [freqVals,psdData,firingRate,dPrimePSD,dPrimeFiringRate,matPSDData,matFiringRate,delPSDData,aucSEPower,aucSEFiringRate] = aucAnalysisPower
folderSourceString = 'G:';
oriType = 'Selected'; % Selected or All
targetOnTimeCriteria = 1; % 0: All trials are considered ; 1:Trials having onset time greater than cutoff is considered
signalDuration = [-0.5 0.1]; % [-1 0.1]
targetOnTimeCutoff = 750; 
folderName = fullfile(folderSourceString,'Projects','MayoProject2','Data','savedDataForAnalysis');
fileNameSaveString = ['data_' oriType 'Ori' '_Perform50' '_targetOnTimeCriterion' num2str(targetOnTimeCriteria) '_targetOnTimeCutoff' num2str(targetOnTimeCutoff) '_targetOnset' num2str(signalDuration(1)) '_' num2str(signalDuration(2)) '.mat']; 
fileName = fullfile(folderName,fileNameSaveString);
data = load(fileName); 

% set up MT parameters
params.tapers = [3 5];
params.pad = -1;
params.Fs = 2000;
params.fpass = [0 200];
params.trialave = 0;

nSessions = 25;
nConditions = 8;
trialCutOff = 10;
timeRange = [-0.5 0];
timePos = intersect(find(data.timeVals>=timeRange(1)),find(data.timeVals<timeRange(2)));

sessionNum = zeros(4,nSessions);
psdData = cell(nConditions,2);
matPSDData = cell(nConditions,2);
delPSDData = cell(4,2);
dPrimePSD = cell(4,2);
aucSEPower = cell(4,2);

firingRate = cell(nConditions,2);
matFiringRate = cell(nConditions,2);
dPrimeFiringRate = cell(4,2);
aucSEFiringRate = cell(4,2);

for c=1:4
    sPrime = 0;
    for s=1:nSessions
        if data.nTrials(s,c)<trialCutOff || data.nTrials(s,c+4)<trialCutOff
            continue
        else
            sPrime = sPrime+1;
            sessionNum(c,sPrime) = s;
            for array=1:2
                nElectrodes = size(data.goodLFPData{array}{s,c},1);
                for e=1:nElectrodes
                    [psdData{c,array}{sPrime}(e,:,:),freqVals] = mtspectrumc(squeeze(data.goodLFPData{array}{s,c}(e,:,timePos))',params);
                    [psdData{c+4,array}{sPrime}(e,:,:)] = mtspectrumc(squeeze(data.goodLFPData{array}{s,c+4}(e,:,timePos))',params);
                    
                    firingRate{c,array}{sPrime}(e,:) = getSpikeCounts(data.goodSpikeData{array}{s,c}(e,:),timeRange)./diff(timeRange);
                    firingRate{c+4,array}{sPrime}(e,:) = getSpikeCounts(data.goodSpikeData{array}{s,c+4}(e,:),timeRange)./diff(timeRange);
                    
                    %ROC analysis
                    clear aucSEPowerTMP
                    for f=1:size(psdData{c,array}{sPrime},2)
                        aucSEPowerTMP(f) = ROCAnalysis(squeeze(psdData{c,array}{sPrime}(e,f,:)),squeeze(psdData{c+4,array}{sPrime}(e,f,:))); %#ok<AGROW>
                    end
                    aucSEPower{c,array} = cat(1,aucSEPower{c,array},aucSEPowerTMP);
                    aucSEFiringRate{c,array} = cat(1,aucSEFiringRate{c,array},ROCAnalysis(firingRate{c,array}{sPrime}(e,:),firingRate{c+4,array}{sPrime}(e,:)));
                end
                dPrimePSD{c,array} = cat(1,dPrimePSD{c,array},abs(mean(psdData{c,array}{sPrime},3) - mean(psdData{c+4,array}{sPrime},3))./sqrt(((data.nTrials(s,c)-1)*var(psdData{c,array}{sPrime},[],3)+(data.nTrials(s,c+4)-1)*var(psdData{c+4,array}{sPrime},[],3))./(data.nTrials(s,c)+data.nTrials(s,c+4)-2)));
                dPrimeFiringRate{c,array} = cat(1,dPrimeFiringRate{c,array},abs(mean(firingRate{c,array}{sPrime}(e,:),2)-mean(firingRate{c+4,array}{sPrime},2))./sqrt(((data.nTrials(s,c)-1)*var(firingRate{c,array}{sPrime},[],2)+(data.nTrials(s,c+4)-1)*var(firingRate{c+4,array}{sPrime},[],2))./(data.nTrials(s,c)+data.nTrials(s,c+4)-2)));
                matPSDData{c,array} = cat(1,matPSDData{c,array},mean(psdData{c,array}{sPrime},3));
                matPSDData{c+4,array} = cat(1,matPSDData{c+4,array},mean(psdData{c+4,array}{sPrime},3));
                matFiringRate{c,array} = cat(1,matFiringRate{c,array},mean(firingRate{c,array}{sPrime},2));
                matFiringRate{c+4,array} = cat(1,matFiringRate{c+4,array},mean(firingRate{c+4,array}{sPrime},2));
                delPSDData{c,array} = cat(1,delPSDData{c,array},log10(mean(psdData{c,array}{sPrime},3))-log10(mean(psdData{c+4,array}{sPrime},3)));
            end
        end
    end
    
end
aucSEPower(:,3) = combineData(aucSEPower);
aucSEFiringRate(:,3) = combineData(aucSEFiringRate);
dPrimePSD(:,3) = combineData(dPrimePSD);
dPrimeFiringRate(:,3) = combineData(dPrimeFiringRate);
matPSDData(:,3) = combineData(matPSDData);
matFiringRate(:,3) = combineData(matFiringRate);
delPSDData(:,3) = combineData(delPSDData);

end

function combinedData = combineData(data)
if size(data,1)==4
    combinedData{1} = cat(1,data{1,1},data{2,2}); %AIV
    combinedData{2} = cat(1,data{1,2},data{2,1}); %AOV
    combinedData{3} = cat(1,data{3,1},data{4,2}); %TIN
    combinedData{4} = cat(1,data{3,2},data{4,1}); %TON
elseif size(data,1)==8
    combinedData{1} = cat(1,data{1,1},data{2,2}); %AIVH
    combinedData{2} = cat(1,data{1,2},data{2,1}); %AOVH
    combinedData{3} = cat(1,data{3,1},data{4,2}); %TIH
    combinedData{4} = cat(1,data{3,2},data{4,1}); %TOH
    combinedData{5} = cat(1,data{5,1},data{6,2}); %AIVM
    combinedData{6} = cat(1,data{5,2},data{6,1}); %AOVM
    combinedData{7} = cat(1,data{7,1},data{8,2}); %TIM
    combinedData{8} = cat(1,data{7,2},data{8,1}); %TOM
end
end

