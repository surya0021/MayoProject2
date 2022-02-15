function [fileNameStringList,monkeyNameList] = getAttentionExperimentDetails

% Arturo
monkeyNameList{1} = 'arturo';
sessionNumberList{1} = [53 55:63 73:75]; % 54 is discarded because of very few H1V and HOI trials %

% Wiggin
monkeyNameList{2} = 'wiggin';
sessionNumberList{2} = [53 55 57 58 60 62 64 65 69 70 72 73];

fileNameStringList = cell(1,2);
for i=1:2
    tmpFileNameList =cell(1,length(sessionNumberList{i}));
    for j=1:length(sessionNumberList{i})
        tmpFileNameList{j} = ['p' monkeyNameList{i}(1) 'cu' num2str(sessionNumberList{i}(j)) 'A001'];
    end
    fileNameStringList{i} = tmpFileNameList;
end

%%%%%%%%%%%%%%%%%%%%%%% Information from the doc file %%%%%%%%%%%%%%%%%%%%%
% i=1;
% fileNameStringList{i} = 'pacu53A001';
% clear stim0 stim1
% stim0.azi=-8; stim0.ele=-5; stim0.dir=40; stim0.rad=3; stim0.sig=1.2; stim0.sf=1;
% stim1.azi=2; stim1.ele=-4.5; stim1.dir=160; stim1.rad=3; stim1.sig=0.8; stim1.sf=1.2;
% stimPropertiesList{i,1} = stim0;
% stimPropertiesList{i,2} = stim1;
% 
% FILE pacu54A001
% stim0: -5, -2; dir=130, rad=3, sig=0.7, sf=1.1
% stim1: 1.5, -5; dir=80, rad=3, sig=0.9, sf=0.6
% 
% FILE pacu55A001
% stim0: -6.5, -4; dir=100, rad=3, sig=0.8, sf=0.8
% stim1: 2, -5.5; dir=90, rad=3, sig=0.9, sf=1.2
% 
% FILE pacu56A001
% stim0: -7.5, -3; dir=100, rad=3, sig=0.9, sf=0.6
% stim1: 2.5, -5.5; dir=0, rad=3, sig=1.1, sf=0.4
% 
% FILE pacu57A001
% stim0: -5, -2.5; dir=0, rad=3, sig=0.9, sf=1.1
% stim1: 2.5, -5; dir=30, rad=3, sig=1.1, sf=1
% 
% FILE pacu58A001
% stim0: -8.5, -3; dir=100, rad=3, sig=0.95, sf=0.6
% stim1: 1.5, -4.5; dir=60, rad=3, sig=0.8, sf=0.5
% 
% FILE pacu59A001
% stim0: -5.5, -2.5; dir=100, rad=3, sig=0.9, sf=0.6
% stim1: 2, -4.5; dir=30, rad=3, sig=0.9, sf=0.4
% 
% FILE pacu60A001
% stim0: -7, -2; dir=90, rad=3, sig=1, sf=0.6
% stim1: 2, -4.5; dir=30, rad=3, sig=1.1, sf=0.9
% 
% FILE pacu61A001
% stim0: -5, -3.5; dir=95, rad=3, sig=0.8, sf=1.2
% stim1: 2, -4; dir=40, rad=3, sig=0.9, sf=1.1
% 
% FILE pacu62A001
% stim0: -6, -4; dir=80, rad=3, sig=1, sf=0.8
% stim1: 2, -7; dir=160, rad=3, sig=1, sf=0.7
% 
% FILE pacu63A001
% stim0: -7, -2.5; dir=130, rad=3, sig=0.8, sf=0.4
% stim1: 2, -5; dir=60, rad=3, sig=1, sf=0.5
% 
% FILE pacu73A001
% stim0: -5.5, -2.5; dir=10, rad=4, sig=0.9, sf=1.2
% stim1: 1.5, -7; dir=0, rad=3, sig=1, sf=0.7
% 
% FILE pacu74A001
% stim0: -5, -1.5; dir=130, rad=4, sig=0.9, sf=0.5
% stim1: 2, -5; dir=30, rad=4, sig=1.1, sf=0.5
% 
% FILE pacu75A001
% stim0: -6, -3; dir=40, rad=4, sig=1, sf=0.6
% stim1: 2, -6; dir=150, rad=4, sig=1, sf=0.7
% 
% 
% 
% 
% 
% 
% FILE pwcu53A001.dat
% stim0: -5.5, -4; dir=122, rad=4, sigma=1, sf=1.7
% stim1: 6.5, -6.5; dir=70, rad=5, sigma=0.7, sf=1.2
% 
% FILE pwcu55A001.dat
% stim0: -5.5, -1.5; dir=30, rad=4, sigma=0.8, sf=0.45
% stim1: 4.5, -7; dir=35, rad=5, sigma=1.1, sf=0.4
% 
% FILE pwcu57A001.dat
% stim0: -7.5, -2.5; dir=78, rad=4, sigma=0.8, sf=1.9
% stim1: 5, -7; dir=45, rad=5, sigma=1.4, sf=2.5
% 
% FILE pwcu58A001.dat
% stim0: -4.5, -1.5; dir=100, rad=4, sigma=0.7, sf=1
% stim1: 4.5, -8.5; dir=90, rad=5, sigma=0.45, sf=2
% 
% FILE pwcu60A001.dat, Cued Uncued.
% 
% FILE pwcu62A001.dat
% stim0: -6.5, -4; dir=130, rad=4, sigma=0.8, sf=0.4
% stim1: 5, -7.5; dir=90, rad=5, sigma=1.5, sf=0.5
% 
% FILE pwcu64A001.dat 
% stim0: -4.5, -2; dir=90, rad=4, sigma=1.1, sf=1.6
% stim1: 5, -7; dir=95, rad=5, sigma=1.6, sf=0.6
% 
% FILE pwcu65A001.dat
% stim0: -7, -2.5; dir=120, rad=4, sigma=1, sf=1.2
% stim1: 6, -8; dir=80, rad=4, sigma=0.6, sf=0.5
% 
% FILE pwcu69A001, Cued Uncued 
% 
% FILE pwcu70A001.dat
% stim0: -7, -4.5; dir=36, rad=4, sigma=1, sf=0.6
% stim1: 3, -9; dir=160, rad=5, sigma=1.6, sf=0.4
% 
% FILE pwcu72A001.dat,
% stim0: -6, -1.5; dir=160, rad=4, sigma=1.1, sf=0.55
% stim1: 3, -6.5; dir=100, rad=5, sigma=1.3, sf=0.3
% 
% FILE pwcu73A001.dat
% stim0: -5, -3.5; dir=85, rad=4, sigma=1, sf=2.2
% stim1: 5, -9; dir=30, rad=4, sigma=0.5, sf=1.6

end