% This program performs ROC analysis given hits and miss condition inputs 

% Improvements required:
% 1)Select equal number of trials in both classes and do a bootstrap without
% replacement
% 2)Check the dependency of number of criterion or distance between
% successive criterion on AUC 


function auc = ROCAnalysis(data1,data2)

criterion = linspace(min(min(data1),min(data2)),max(max(data1),max(data2)),100);

class1 = data1;
class2 = data2;

for i=1:length(criterion)
    TP(i) = sum(class1>=criterion(i)); %#ok<*AGROW>
    FN(i) = sum(class1<criterion(i));
    TN(i) = sum(class2<criterion(i));
    FP(i) = sum(class2>=criterion(i));
    
    TPR(i) = TP(i)/(TP(i)+FN(i));
    FPR(i) = FP(i)/(FP(i)+TN(i));
end
auc = abs(trapz(FPR,TPR));
end