%% Universal Settings for Building Populations 

% cardiomyocyte stimulation settings 
settings.celltype = 'endo';
settings.PCL = 1000 ;  % Pacing, Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 5 ; % Stimulus duration
settings.stim_amp = 12.5; % Stimulus amplitude 

% variability settings 
settings.variations = 100000; % population size 
settings.sigmaG = 0.4; % standard deviation to vary conductances 
% calibrate population? 
% when creating the initial population, sometimes certain parameter sets
% create cells that form arrhythmic activity before any trigger is applied.
% we removed those cells, and reran new ones.
settings.remove_arrhythmias = true; % remove cells that form arrhythmic activity 
settings.remove_experimental = false; % remove cells not within experimental range
settings.reruncells = true; % rerun removed cells 

%-- Folder to Save Data 
settings.Folder = 'TestPop';
if ~exist(settings.Folder, 'dir')
    mkdir(settings.Folder)
end
addpath('TestPop');
%-- Run Population 
pert = settings_blockcurrents; 
BuildPopulation(settings,pert)
disp('Saving population  to TestPop/BaseCells')