clear
load('TestPop\SR\male\SRpopm.mat') % SR male population data
load('TestPop\SR\female\SRpopf.mat') % SR female population data

%% Steady State Features
% SR male
for j = 1:length(SRBaseCellsm)
    t = SRBaseCellsm(j).times;
    V = SRBaseCellsm(j).V;
    Cai = SRBaseCellsm(j).Cai .* 1000000;
    SRfeaturesm(j,:) = calculate_features(V,Cai,t);    
end
SRfeaturesm(:,3) = SRfeaturesm(:,3)-20;

% SR female
for j = 1:length(SRBaseCellsf)
    t = SRBaseCellsf(j).times;
    V = SRBaseCellsf(j).V;
    Cai = SRBaseCellsf(j).Cai .* 1000000;
    SRfeaturesf(j,:) = calculate_features(V,Cai,t);    
end
SRfeaturesf(:,3) = SRfeaturesf(:,3)-15;

fetureName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'CTA', 'CTmax', 'CTD50', 'CTD90', 'CTDtri', 'CTD'};

%% Save feature date 
matfile = fullfile('TestPop\SR\male\', 'SRfeaturesm.mat');
save(matfile,'SRfeaturesm','fetureName')

matfile = fullfile('TestPop\SR\female\', 'SRfeaturesf.mat');
save(matfile,'SRfeaturesf','fetureName')

disp('Saving population feature')
%}