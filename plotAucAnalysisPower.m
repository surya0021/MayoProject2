% This program plots PSD, firing rate, dPrime anf AUC of single electrodes and averaged across electrodes for each attention
% condition seperately.

% get data
[freqVals,psdData,firingRate,dPrimePSD,dPrimeFiringRate,matPSDData,matFiringRate,delPSDData,aucSEPower,aucSEFiringRate] = aucAnalysisPower;
conditionList = {'Attend-in Valid','Attend-out Valid','Target-in Neutral','Target-out Neutral'};
logFlag=1;
for c=1:4
    %Figure
    set(matlab.ui.Figure,'units','normalized','outerposition',[0 0 1 1]);
    hSEMeasures1 = getPlotHandles(1,1,[0.05 0.65 0.1 0.25],0.1,0.1,0);
    hSEMeasures = getPlotHandles(1,3,[0.2 0.65 0.75 0.25],0.05,0.05,0);
    hAvgMeasures1 = getPlotHandles(1,1,[0.05 0.35 0.1 0.25],0.1,0.1,0);
    hAvgMeasures = getPlotHandles(1,3,[0.2 0.35 0.75 0.25],0.05,0.05,0);
    hAvgMeasures2 = getPlotHandles(1,1,[0.3 0.05 0.5 0.25],0.05,0.05,0);
    sgtitle(conditionList{c})
    % Plotting single electrode measures
    cc = c; ac = 1; sc = 1; ec = 1; ee = 1;
    % Firing rate
    C = categorical({'Hits','Miss'});
    matFR = [mean(firingRate{cc,ac}{sc}(ec,:),2),mean(firingRate{cc+4,ac}{sc}(ec,:),2)];
    sematFR = [std(firingRate{cc,ac}{sc}(ec,:),[],2)/sqrt(size(firingRate{cc,ac}{sc},2)),std(firingRate{cc+4,ac}{sc}(ec,:),[],2)/sqrt(size(firingRate{cc+4,ac}{sc},2))];
    a=bar(hSEMeasures1,C,matFR,0.5);
    hold(hSEMeasures1,'on')
    errorbar(hSEMeasures1,[1,2],matFR,sematFR,'color','k','lineStyle','none');
    ylabel(hSEMeasures1,'Firing rate (spikes/s)')
    text(hSEMeasures1,a.XData(2),45,['dPrime:' num2str(dPrimeFiringRate{1,1}(33))])
    
    % PSD
    if logFlag==1
        matPSD = {mean(log10(psdData{cc,ac}{sc}(ec,:,:)),3),mean(log10(psdData{cc+4,ac}{sc}(ec,:,:)),3)};
        sematPSD = {std(log10(psdData{cc,ac}{sc}(ec,:,:)),[],3)./sqrt(size(psdData{cc,ac}{sc},3)),std(log10(psdData{cc+4,ac}{sc}(ec,:,:)),[],3)./sqrt(size(psdData{cc+4,ac}{sc},3))};
    elseif logFlag==0
        matPSD = {mean((psdData{cc,ac}{sc}(ec,:,:)),3),mean((psdData{cc+4,ac}{sc}(ec,:,:)),3)};
        sematPSD = {std((psdData{cc,ac}{sc}(ec,:,:)),[],3)./sqrt(size(psdData{cc,ac}{sc},3)),std((psdData{cc+4,ac}{sc}(ec,:,:)),[],3)./sqrt(size(psdData{cc+4,ac}{sc},3))};
    end
     % Hit
    p1 = plot(hSEMeasures(1,1),freqVals,matPSD{1},'b');
    hold(hSEMeasures(1,1),'on')
    patch([freqVals fliplr(freqVals)],[matPSD{1}+sematPSD{1} fliplr(matPSD{1}-sematPSD{1})],'b','FaceAlpha',0.5,'lineWidth',1,'parent',hSEMeasures(1,1))
    %Miss
    p2 = plot(hSEMeasures(1,1),freqVals,matPSD{2},'r');
    patch([freqVals fliplr(freqVals)],[matPSD{2}+sematPSD{2} fliplr(matPSD{2}-sematPSD{2})],'r','FaceAlpha',0.5,'lineWidth',1,'parent',hSEMeasures(1,1))
    xlabel(hSEMeasures(1,1),'Frequency (Hz)')
    ylabel(hSEMeasures(1,1), 'Log_{10} power (\muV)');
    legend(hSEMeasures(1,1),[p1 p2],{'Hit','Miss'},'location','northeast');     legend(hSEMeasures(2),'boxoff')
    
    % dPrime and AUC
    plot(hSEMeasures(1,2),freqVals,dPrimePSD{cc,ac}(ee,:),'b')
    hold(hSEMeasures(1,2),'on')
    plot(hSEMeasures(1,2),freqVals,repmat(dPrimeFiringRate{cc,ac}(ee),1,length(freqVals)),'r')
    xlabel(hSEMeasures(1,2),'Frequency (Hz)')
    ylabel(hSEMeasures(1,2),'d-Prime')
    
    plot(hSEMeasures(1,3),freqVals,aucSEPower{cc,ac}(ee,:),'b')
    hold (hSEMeasures(1,3),'on')
    plot(hSEMeasures(1,3),freqVals,repmat(aucSEFiringRate{cc,ac}(ee,:),1,length(freqVals)),'r')
    xlabel(hSEMeasures(1,3),'Frequency (Hz)')
    ylabel(hSEMeasures(1,3),'AUC')
    
    
    % For average measures
    
    % Firing rate
    C = categorical({'Hits','Miss'});
    maeFR = [mean(matFiringRate{cc,3},1),mean(matFiringRate{cc+4,3},1)];
    semaeFR = [std(matFiringRate{cc,3},[],1)/sqrt(size(matFiringRate{cc,3},1)),std(matFiringRate{cc+4,3},[],1)/sqrt(size(matFiringRate{cc+4,3},1))];
    a=bar(hAvgMeasures1,C,maeFR,0.5);
    hold(hAvgMeasures1,'on')
    errorbar(hAvgMeasures1,[1,2],maeFR,semaeFR,'color','k','lineStyle','none');
    ylabel(hAvgMeasures1,'Firing rate (spikes/s)')
    
    
    
    % PSD
    if logFlag==1
        maePSD = {mean(log10(matPSDData{cc,3}),1),mean(log10(matPSDData{cc+4,3}),1)};
        semaePSD = {std(log10(matPSDData{cc,3}),[],1)./sqrt(size(matPSDData{cc,3},1)),std(log10(matPSDData{cc+4,3}),[],1)./sqrt(size(matPSDData{cc+4,3},1))};
    elseif logFlag==0
        maePSD = {mean((matPSDData{cc,3}),1),mean((matPSDData{cc+4,3}),1)};
        semaePSD = {std((matPSDData{cc,3}),[],1)./sqrt(size(matPSDData{cc,3},1)),std((matPSDData{cc+4,3}),[],1)./sqrt(size(matPSDData{cc+4,3},1))};
    end
    % Hit
    p1 = plot(hAvgMeasures(1,1),freqVals,maePSD{1},'b');
    hold(hAvgMeasures(1,1),'on')
    patch([freqVals fliplr(freqVals)],[maePSD{1}+semaePSD{1} fliplr(maePSD{1}-semaePSD{1})],'b','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures(1,1))
    %Miss
    p2 = plot(hAvgMeasures(1,1),freqVals,maePSD{2},'r');
    patch([freqVals fliplr(freqVals)],[maePSD{2}+semaePSD{2} fliplr(maePSD{2}-semaePSD{2})],'r','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures(1,1))
    xlabel(hAvgMeasures(1,1),'Frequency (Hz)')
    ylabel(hAvgMeasures(1,1), 'Log_{10} power (\muV)');
    legend(hAvgMeasures(1,1),[p1 p2],{'Hit','Miss'},'location','northeast');     legend(hAvgMeasures(2),'boxoff')
    
    
    
    % dPrime
    maeDPrimePSD = mean(dPrimePSD{cc,3},1);
    semaeDPrimePSD = std(dPrimePSD{cc,3},[],1)./sqrt(size(dPrimePSD{cc,3},1));
    maeDPrimeFiringRate = mean(dPrimeFiringRate{cc,3},1);
    semaeDPrimeFiringRate = std(dPrimeFiringRate{cc,3},[],1)/sqrt(size(dPrimeFiringRate{cc,3},1));
    
    plot(hAvgMeasures(1,3),freqVals,maeDPrimePSD,'b')
    hold(hAvgMeasures(1,3),'on')
    patch([freqVals fliplr(freqVals)],[maeDPrimePSD+semaeDPrimePSD fliplr(maeDPrimePSD-semaeDPrimePSD)],'b','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures(1,3))
    plot(hAvgMeasures(1,3),freqVals,repmat(maeDPrimeFiringRate,1,length(freqVals)),'r')
    patch([freqVals fliplr(freqVals)],[repmat(maeDPrimeFiringRate+semaeDPrimeFiringRate,1,length(freqVals)) repmat(maeDPrimeFiringRate-semaeDPrimeFiringRate,1,length(freqVals))],'r','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures(1,3));
    xlabel(hAvgMeasures(1,3),'Frequency (Hz)')
    ylabel(hAvgMeasures(1,3),'d-Prime')
    
    
    
    % Change in Power
    if logFlag==1
        maeDelPSDData = mean(delPSDData{cc,3},1);
        semaeDelPSDData = std(delPSDData{cc,3},[],1)./sqrt(size(delPSDData{cc,3},1));
    elseif logFlag==0
        maeDelPSDData = mean(matPSDData{cc,3}-matPSDData{cc+4,3},1);
        semaeDelPSDData = std(matPSDData{cc,3}-matPSDData{cc+4,3},[],1)./sqrt(size(matPSDData{cc,3},1));
    end
    plot(hAvgMeasures(1,2),freqVals,maeDelPSDData,'b')
    hold(hAvgMeasures(1,2),'on')
    plot(hAvgMeasures(1,2),freqVals,zeros(1,length(freqVals)),'k')
    patch([freqVals fliplr(freqVals)],[maeDelPSDData+semaeDelPSDData fliplr(maeDelPSDData-semaeDelPSDData)],'b','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures(1,2))
    xlabel(hAvgMeasures(1,2),'Frequency (Hz)')
    ylabel(hAvgMeasures(1,2),'Change in power')
    % AUC
    maeAucSEPower = mean(aucSEPower{cc,3},1);
    semaeAucSEPower = std(aucSEPower{cc,3},[],1)./sqrt(size(aucSEPower{cc,3},1));
    maeAucSEFiringRate = mean(aucSEFiringRate{cc,3},1);
    semaeAucSEFiringRate = std(aucSEFiringRate{cc,3},[],1)/sqrt(size(aucSEFiringRate{cc,3},1));
    
    plot(hAvgMeasures2,freqVals,maeAucSEPower,'b')
    hold (hAvgMeasures2,'on')
    patch([freqVals fliplr(freqVals)],[maeAucSEPower+semaeAucSEPower fliplr(maeAucSEPower-semaeAucSEPower)],'b','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures2);
    plot(hAvgMeasures2,freqVals,repmat(maeAucSEFiringRate,1,length(freqVals)),'r')
    patch([freqVals fliplr(freqVals)],[repmat(maeAucSEFiringRate+semaeAucSEFiringRate,1,length(freqVals)) repmat(maeAucSEFiringRate-semaeAucSEFiringRate,1,length(freqVals))],'r','FaceAlpha',0.5,'lineWidth',1,'parent',hAvgMeasures2);
    xlabel(hAvgMeasures2,'Frequency (Hz)')
    ylabel(hAvgMeasures2,'AUC')
    
end