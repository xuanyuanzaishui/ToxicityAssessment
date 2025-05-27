function BuildPopulation(settings,pert)
%--------------------------------------------------------------------------
%% -- BuildPopulation.m -- 
% Description: Build population

% Inputs:
% --> settings - [struct] simulation protocol settings
% --> pert  - [struct] initialize channel block settings

% Outputs: No actual outputs. Data is automatically saved to folder.
% -------------------------------------------------------------------------
%% Settings
settings.nBeats = 500; % Number of beats to simulate
settings.numbertokeep = 2 ;% Determine how many beats to keep. 2 = last two beat
settings.steady_state = true; % Start with each scenario's BL model steady state values.

%% Run Simulation

disp('Calibrate Population...')

% Does the user provide a parameter matrix or should we make one?
if ~isfield(settings,'scalings') % 检查settings是否包含scalings，包含0，不包含1
    settings.scalings = popfuncs.create_scale_vector(settings,settings.variations);
end

% Separate population
intervals = chunks(settings,settings.variations);
chunk_settings = settings;
% define var
popICs = [];
popscalings = [];
AFpopICs = [];
AFpopscalings = [];
a=1;
b=1;
% Loop through multiple intervals of the data
for ii = 1:size(intervals,1)
    chunk_settings.scalings = settings.scalings(intervals(ii,1):intervals(ii,2),:);
    [chunk_settings.variations,~] = size(chunk_settings.scalings);

    X = runSim(chunk_settings,pert); % run simulation
    Y = popfuncs.clean_population(X,chunk_settings);% clean data

    if ~isempty(Y)
      for i = 1:length(Y)
        BaseCells(a).times = Y(i).times;
        BaseCells(a).V = Y(i).states(:,39);
        BaseCells(a).Cai = Y(i).states(:,38);
        popICs(end+1,:) = Y(i).states(end,:);
        popscalings(end+1,:) = Y(i).scalings;
        a=a+1;
      end
    end

    clear X Y
    disp([num2str((ii/(size(intervals,1)))*100) '% Finished '])
end

matfile = fullfile(settings.Folder, 'BasePop.mat');
save(matfile, 'BaseCells','popICs','popscalings')

%% -- Nested Functions 
function intervals = chunks(~,variations)
    if variations > 1000
        ints = round(linspace(0,variations,variations/500));
        intervals = zeros((length(ints) - 1),2);
        for i = 1:(length(ints) - 1)
            intervals(i,1:2) = [ints(i)+1 ints(i+1)];
        end
    else
        intervals = [1 variations];
    end