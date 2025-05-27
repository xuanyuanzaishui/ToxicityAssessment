clear

%% Universal Settings for Calibrate Populations 
load("basepop.mat")
load("basepop_init.mat")

BaseCell(1).times = BaseCells(1).times;
BaseCell(1).V = BaseCells(1).V;
BaseCell(1).Cai = BaseCells(1).Cai;

c=1;
d=1;
for i = 1:length(BaseCells)
    t = BaseCells(i).times;
    V = BaseCells(i).V;
    Cai = BaseCells(i).Cai .* 1000000;
    a=1;
    for ii = 1:length(t)
        if V(ii) < (max(V)-21.5)
           ti(a,1) = t(ii); 
           Vi(a,1) = V(ii); 
           Caii(a,1) = Cai(ii); 
           a=a+1;
        end
    end

    outputs = calculate_features(Vi,Caii,ti);

    RMP = outputs(1,1);
    dVdt = outputs(1,2);
    APA = outputs(1,3);
    APD20 = outputs(1,4);
    APD50 = outputs(1,6);
    APD90 = outputs(1,7);
    
    %AFmale
    I1 = RMP > (-82.12) && RMP < (-73.02);
    I2 = APA > 93.6 && APA < 113.8;
    I3 = dVdt > 177.26 && dVdt < 339.2;
    I4 = APD20 > 11.23 && APD20 < 50;
    I5 = APD50 > 68.29 && APD50 < 130.09;
    I6 = APD90 > 172.21 && APD90 < 269.81;
    I = [I1,I2,I3,I4,I5,I6];
    %AFfemale
    J1 = RMP > (-81.82) && RMP < (-72.72);
    J2 = APA > 92.38 && APA < 112.58;
    J3 = dVdt > 162.13 && dVdt < 324.07;
    J4 = APD20 > 12.2 && APD20 < 51;
    J5 = APD50 > 70.89 && APD50 < 132.69;
    J6 = APD90 > 178.39 && APD90 < 275.99;   
    J = [J1,J2,J3,J4,J5,J6];

    if all(I == 1)
        AFBaseCellsm(c).times = BaseCells(i).times;
        AFBaseCellsm(c).V = BaseCells(i).V;
        AFBaseCellsm(c).Cai = BaseCells(i).Cai;
        AFpopICsm(c,:) = popICs(i,:);
        AFpopscalingsm(c,:) = popscalings(i,:);
        c=c+1;
    end
    if all(J == 1)
        AFBaseCellsf(d).times = BaseCells(i).times;
        AFBaseCellsf(d).V = BaseCells(i).V;
        AFBaseCellsf(d).Cai = BaseCells(i).Cai;
        AFpopICsf(d,:) = popICs(i,:);
        AFpopscalingsf(d,:) = popscalings(i,:);
        d=d+1;
    end
           
end

matfile = fullfile('TestPop\AF\male\', 'AFpopICsm.mat');
save(matfile,'AFpopICsm','AFpopscalingsm')
matfile = fullfile('TestPop\AF\male\', 'AFpopm.mat');
save(matfile,'AFBaseCellsm')

matfile = fullfile('TestPop\AF\female\', 'AFpopICsf.mat');
save(matfile,'AFpopICsf','AFpopscalingsf')
matfile = fullfile('TestPop\AF\female\', 'AFpopf.mat');
save(matfile,'AFBaseCellsf')

disp('Saving population')
