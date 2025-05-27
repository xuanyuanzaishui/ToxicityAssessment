clear
load('TestPop\SR2014\SRpop.mat') % SR male population data
load('TestPop\SR2014\SRpopICs.mat') % SR male population data


%% Steady State Features
% SR male
a=1;
for j = 1:length(SRBaseCells2014)
    t = SRBaseCells2014(j).times;
    V = SRBaseCells2014(j).V;
    Cai = SRBaseCells2014(j).Cai .* 1000000;
    SRfeatures(j,:) = calculate_features(V,Cai,t); 

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
        SRBaseCells20141(a).times = SRBaseCells2014(j).times;
        SRBaseCells20141(a).V = SRBaseCells2014(j).V;
        SRBaseCells20141(a).Cai = SRBaseCells2014(j).Cai;
        SRpopICs20141(a,:) = SRpopICs2014(j,:);
        SRpopscalings20141(a,:) = SRpopscalings2014(j,:);
        SRfeatures2014(a,:) = SRfeatures(j,:);
        a=a+1;
    end
end
SRfeatures2014(:,3) = SRfeatures2014(:,3)-10;

mean1 = mean(SRfeatures2014);
std1 = std(SRfeatures2014);
meanstd = [mean1;std1];

fetureName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'CTA', 'CTmax', 'CTD50', 'CTD90', 'CTDtri', 'CTD'};

%% Save feature date 
matfile = fullfile('TestPop\SR2014\', 'SRfeatures2014.mat');
save(matfile,'SRfeatures2014','fetureName')
matfile = fullfile('TestPop\SR2014\', 'SRpopICs2014.mat');
save(matfile,'SRpopICs20141','SRpopscalings20141') 
matfile = fullfile('TestPop\SR2014\', 'SRpop2014.mat');
save(matfile,'SRBaseCells20141')

disp('Saving population feature')
%}