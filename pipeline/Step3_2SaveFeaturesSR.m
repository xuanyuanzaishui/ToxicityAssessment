clear
load('TestPop\SR2022\SRpop.mat') % SR male population data
load('TestPop\SR2022\SRpopICs.mat') % SR male population data


%% Steady State Features
% SR male
a=1;
for j = 1:length(SRBaseCells2022)
    t = SRBaseCells2022(j).times;
    V = SRBaseCells2022(j).V;
    Cai = SRBaseCells2022(j).Cai .* 1000000;
    SRfeatures(j,:) = calculate_features(V,Cai,t); 
%{
    RMP = SRfeatures(j,1);
    dVdt = SRfeatures(j,2);
    APA = SRfeatures(j,3)-10;
    APD20 = SRfeatures(j,4);
    APD50 = SRfeatures(j,6);
    APD90 = SRfeatures(j,7);
    %SR2014
    II1 = RMP > (-85) && RMP < (-65);
    II2 = APA > 75 && APA < 120;
    II3 = dVdt > 135 && dVdt < 370;
    II4 = APD20 > 1 && APD20 < 60;
    II5 = APD50 > 75 && APD50 < 200;
    II6 = APD90 > 255 && APD90 < 420;
    II = [II1,II2,II3,II4,II5,II6];
    if all(II == 1)
        SRBaseCells20221(a).times = SRBaseCells2022(j).times;
        SRBaseCells20221(a).V = SRBaseCells2022(j).V;
        SRBaseCells20221(a).Cai = SRBaseCells2022(j).Cai;
        SRpopICs20221(a,:) = SRpopICs2022(j,:);
        SRpopscalings20221(a,:) = SRpopscalings2022(j,:);
        SRfeatures2022(a,:) = SRfeatures(j,:);
        a=a+1;
    end
%}
end
SRfeatures2022 = SRfeatures;
SRfeatures2022(:,3) = SRfeatures2022(:,3)-10;

mean1 = mean(SRfeatures2022);
std1 = std(SRfeatures2022);
meanstd = [mean1;std1];

fetureName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'CTA', 'CTmax', 'CTD50', 'CTD90', 'CTDtri', 'CTD'};

%% Save feature date 
matfile = fullfile('TestPop\SR2022\', 'SRfeatures2022.mat');
save(matfile,'SRfeatures2022','fetureName')

disp('Saving population feature')
%}