function ApplyTrigger(settings,pert)

%% Standard Settings
settings.nBeats = 100 ; % Number of beats to simulate
settings.numbertokeep = 2 ;% Determine how many beats to keep. 1 = last beat, 2 = last two beats

%% Load ICs from Folder 
% SR population data
%File = fullfile('TestPop\SR\male\SRpopICsm.mat'); % male 
File = fullfile('TestPop\SR\female\SRpopICsf.mat'); % female 


% AF population data
%File = fullfile('TestPop\AF\male\AFpopICsm.mat'); % male 
%File = fullfile('TestPop\AF\female\AFpopICsf.mat'); % female

%Modify loading parameter values based on loaded population data
load(File,'SRpopscalings','SRpopICs');
[variations,~] = size(SRpopscalings);


cx=1;%用来矩阵计数
CApopICs = [];
CApopscalings = [];

%% Set Intervals 
% Since populations are so large, we need to divide how we save the data
% into intervals. Number of sections is determined by input.
intervals = chunks(variations);
chunk_settings = settings;
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
if ~exist(yourFolder, 'dir')
    mkdir(yourFolder)
end

%% Run Simulation
for ii = 1:size(intervals,1)
    % Set Up Population Variants
    n = intervals(ii,1):intervals(ii,2);
    chunk_settings.ICs = SRpopICs(n,:);
    chunk_settings.scalings = SRpopscalings(n,:);
    chunk_settings.variations = length(n);
    
    X = runSim(chunk_settings,pert); % run simulation    
    
    for jj = 1:length(X) %for each member of the population
        CABaseCells(cx).times = X(jj).times;
        CABaseCells(cx).V = X(jj).states(:,39);
        CABaseCells(cx).Cai = X(jj).states(:,38);
        CApopICs(end+1,:) = X(jj).states(end,:);
        CApopscalings(end+1,:) = X(jj).scalings;
        cx = cx + 1;
    end
    clear X
    disp([num2str((ii/size(intervals,1))*100) '% Finished '])
end

% Save Data
matfile = fullfile(yourFolder, 'CABaseCells.mat');
save(matfile,'CABaseCells','CApopICs','CApopscalings')

clear CABaseCells

%% -- Nested Functions 
function intervals = chunks(variations)
    if variations > 1000
        ints = round(linspace(0,variations,variations/100));
        intervals = zeros((length(ints) - 1),2);
        for i = 1:(length(ints) - 1)
            intervals(i,1:2) = [ints(i)+1 ints(i+1)];
        end
    else
        intervals = [1 variations];
    end