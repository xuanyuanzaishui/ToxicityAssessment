%% Universal Settings for Applying ClassIII drugs 
clear
% cardiomyocyte settings
settings.celltype = 'endo'; % 'epi', 'endo', 'mid',
settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 5 ; % Stimulus duration
settings.stim_amp = 12.5; % Stimulus amplitude 

%settings.Folder = 'TestPop/SRClassIII/male';
settings.Folder = 'TestPop/SRClassIII/female';
%settings.Folder = 'TestPop/SRClassIII2022';


%settings.Folder = 'TestPop/AFClassIII/male';
%settings.Folder = 'TestPop/AFClassIII/female';

%-- Apply Amiodarone Block Trigger
%
settings.SubFolder = 'Amiodarone';
pert = settings_blockcurrents; 
pert.GNa = SimulateDrug(2,5,1);
pert.Gto = SimulateDrug(2,3.8,0.4);
pert.GKr = SimulateDrug(2,3.0,1);
pert.GKs = SimulateDrug(2,100,1);
pert.GNCX  = SimulateDrug(2,3.4,1);
pert.GCaL = SimulateDrug(2,1.5,0.6);
ApplyTrigger(settings,pert)
pause(25)


%-- Apply Dofetilide Trigger
settings.SubFolder = 'Dofetilide';
pert = settings_blockcurrents; 
pert.GNa = SimulateDrug(2,380.5,0.9);
pert.GKr = SimulateDrug(2,4.9,0.9);
pert.GCaL = SimulateDrug(2,260.3,1.2);
pert.GNaL = SimulateDrug(2,753160.4,0.3);
pert.Gto = SimulateDrug(2,18.8,0.8);
pert.GK1 = SimulateDrug(2,394.3,0.8);

ApplyTrigger(settings,pert)
pause(25)

%-- Apply 
settings.SubFolder = 'Dronedarone';
pert = settings_blockcurrents; 
pert.GNa = SimulateDrug(0.205,0.54,2.03);
pert.GKr = SimulateDrug(0.205,0.0591,0.8);
pert.GKs = SimulateDrug(0.205,5.6,0.51);
pert.GKur = SimulateDrug(0.205,1,1);
pert.GCaL = SimulateDrug(0.205,0.83,2.75);
ApplyTrigger(settings,pert)
pause(25)

%-- Apply Ibutilide Trigger
settings.SubFolder = 'Ibutilide';
pert = settings_blockcurrents; 

pert.GKr = SimulateDrug(0.015,0.02,1);
ApplyTrigger(settings,pert)
pause(25)


%-- Apply Sotalol Trigger
settings.SubFolder = 'Sotalol';
pert = settings_blockcurrents; 
pert.GNa = SimulateDrug(14690,1140000000,0.5);
pert.Gto = SimulateDrug(14690,43143455,0.7);
pert.GKr = SimulateDrug(14690,110600,0.8);
pert.GKs = SimulateDrug(14690,4221856,1.2);
pert.GCaL = SimulateDrug(14690,7061527,0.9);
pert.GK1 = SimulateDrug(14690,3050260,1.2);
ApplyTrigger(settings,pert)
pause(25)
%}

%-- Apply Vernakalant Trigger
settings.SubFolder = 'Vernakalant';
pert = settings_blockcurrents; 
pert.GNa = SimulateDrug(10,90,1);
pert.Gto = SimulateDrug(10,15,1);
pert.GKr = SimulateDrug(10,20,1);
pert.GKur = SimulateDrug(10,15,1);
pert.GCaL = SimulateDrug(10,84,1);
ApplyTrigger(settings,pert)


disp('finished')

%% -- Drugs model Functions 
function Gi_drug = SimulateDrug(Drug,IC50i,hill)
    Gi_drug = 1/(1+((Drug)/IC50i)^(hill));
end