clear
%% Universal Settings for Calibrate Populations 
load("SRBaseCells.mat")
load("SRBaseICs.mat")

BaseCell(1).times = BaseCells(1).times;
BaseCell(1).V = BaseCells(1).V;
BaseCell(1).Cai = BaseCells(1).Cai;

a=1;
b=1;
for i = 1:length(BaseCells)
    ti = BaseCells(i).times;
    Vi = BaseCells(i).V;
    Caii = BaseCells(i).Cai .* 1000000;
    outputs = calculate_features(Vi,Caii,ti);

    RMP = outputs(1,1);
    dVdt = outputs(1,2);
    APA = outputs(1,3);
    APD20 = outputs(1,4);
    APD50 = outputs(1,6);
    APD90 = outputs(1,7);

    APD40 = outputs(1,5);
    CTA = outputs(1,9);
    CTD90 = outputs(1,12);
    %SR2022
    II1 = RMP > (-84) && RMP < (-69);
    II2 = APA > 97 && APA < 135;
    II3 = dVdt > 110 && dVdt < 370;
    II4 = APD20 > 0.4 && APD20 < 15;
    II5 = APD50 > 103 && APD50 < 175;
    II6 = APD90 > 250 && APD90 < 320;
    II = [II1,II2,II3,II4,II5,II6];
    if all(II == 1)
        SRBaseCells2022(a).times = BaseCells(i).times;
        SRBaseCells2022(a).V = BaseCells(i).V;
        SRBaseCells2022(a).Cai = BaseCells(i).Cai;
        SRpopICs2022(a,:) = popICs(i,:);
        SRpopscalings2022(a,:) = popscalings(i,:);
        a=a+1;
    end

    %{
    %SR2014
    II1 = RMP > (-85) && RMP < (-65);
    II2 = APA > 75 && APA < 120;
    II3 = dVdt > 115 && dVdt < 370;
    II4 = APD20 > 1 && APD20 < 60;
    II5 = APD50 > 75 && APD50 < 200;
    II6 = APD90 > 210 && APD90 < 420;
    II = [II1,II2,II3,II4,II5,II6];
    if all(II == 1)
        SRBaseCells2014(a).times = BaseCells(i).times;
        SRBaseCells2014(a).V = BaseCells(i).V;
        SRBaseCells2014(a).Cai = BaseCells(i).Cai;
        SRpopICs2014(a,:) = popICs(i,:);
        SRpopscalings2014(a,:) = popscalings(i,:);
        a=a+1;
    end
    %}
    %{
    %SR2023male
    II1 = RMP > (-77.45) && RMP < (-69.71);
    II2 = APA > 85.3 && APA < 102.7;
    II3 = dVdt > 188.04 && dVdt < 320.04;
    II4 = APD20 > 0 && APD20 < 40;
    II5 = APD50 > 96.08 && APD50 < 187.88;
    II6 = APD90 > 278.11 && APD90 < 385;
    %II7 = APD40 > 70 && APD40 < 114;
    %II8 = CTD90 > 490 && CTD90 < 630;
    %II9 = CTA > 210 && CTA < 420;
    II = [II1,II2,II3,II4,II5,II6];

    %SR2023female
    JJ1 = RMP > (-76.03) && RMP < (-67.3);
    JJ2 = APA > 82.6 && APA < 100;
    JJ3 = dVdt > 176.07 && dVdt < 298.07;
    JJ4 = APD20 > 3 && APD20 < 45;
    JJ5 = APD50 > 105.97 && APD50 < 207.77;
    JJ6 = APD90 > 251 && APD90 < 356.6;
    %JJ7 = APD40> 55 && APD40 < 110;
    %JJ8 = CTD90 > 490 && CTD90 < 650;
    %JJ9 = CTA > 190 && CTA < 400;
    JJ = [JJ1,JJ2,JJ3,JJ4,JJ5,JJ6];
    if all(II == 1)
        SRBaseCellsm(a).times = BaseCells(i).times;
        SRBaseCellsm(a).V = BaseCells(i).V;
        SRBaseCellsm(a).Cai = BaseCells(i).Cai;
        SRpopICsm(a,:) = popICs(i,:);
        SRpopscalingsm(a,:) = popscalings(i,:);
        a=a+1;
    end
    if all(JJ == 1)
        SRBaseCellsf(b).times = BaseCells(i).times;
        SRBaseCellsf(b).V = BaseCells(i).V;
        SRBaseCellsf(b).Cai = BaseCells(i).Cai;
        SRpopICsf(b,:) = popICs(i,:);
        SRpopscalingsf(b,:) = popscalings(i,:);
        b=b+1;
    end
    %}    
end


matfile = fullfile('TestPop\SR2022\', 'SRpopICs.mat');
save(matfile,'SRpopICs2022','SRpopscalings2022')
matfile = fullfile('TestPop\SR2022\', 'SRpop.mat');
save(matfile,'SRBaseCells2022')

%{
matfile = fullfile('TestPop\SR\female\', 'SRpopICsf.mat');
save(matfile,'SRpopICsf','SRpopscalingsf')
matfile = fullfile('TestPop\SR\female\', 'SRpopf.mat');
save(matfile,'SRBaseCellsf')
%}
disp('Saving population')
