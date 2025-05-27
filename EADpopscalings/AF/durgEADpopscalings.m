clear

%% 1、查找EAD  2、保存ICs、scaling


%
%查找EAD male
load('TestPop/AF/male/AFpopICsm.mat','AFpopscalingsm');
load('TestPop/AF/female/AFpopICsf.mat','AFpopscalingsf');

settings.Folder = 'TestPop/AFClassIII/male';%AF
settings.Folder1 = 'TestPop/AFClassIII/female';%AF

%% Dofetilide 
%
settings.SubFolder = 'Dofetilide';
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = AFpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 1;
popscalingsm(1:end,18) = 0;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = AFpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 1;
popscalingsf(1:end,18) = 0;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model21\conduction\'; % path
Output_File1 = 'AFEADpopscalingsmf2.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}


%% Dronedarone
%
settings.SubFolder = 'Dronedarone';

%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = AFpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 1;
popscalingsm(1:end,18) = 0;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = AFpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 1;
popscalingsf(1:end,18) = 0;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model21\conduction\'; % path
Output_File1 = 'AFEADpopscalingsmf3.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}


%% Ibutilide
%
settings.SubFolder = 'Ibutilide';
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = AFpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 1;
popscalingsm(1:end,18) = 0;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = AFpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 1;
popscalingsf(1:end,18) = 0;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model21\conduction\'; % path
Output_File1 = 'AFEADpopscalingsmf4.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}

%% Sotalol
%
settings.SubFolder = 'Sotalol';
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = AFpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 1;
popscalingsm(1:end,18) = 0;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = AFpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 1;
popscalingsf(1:end,18) = 0;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model21\conduction\'; % path
Output_File1 = 'AFEADpopscalingsmf5.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}

%% Vernakalant
%
settings.SubFolder = 'Vernakalant';
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = AFpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 1;
popscalingsm(1:end,18) = 0;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = AFpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 1;
popscalingsf(1:end,18) = 0;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model21\conduction\'; % path
Output_File1 = 'AFEADpopscalingsmf6.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}


%% function
function outputFile(Output_Dir, Output_File, outputFeature)
    % 构建完整路径
    path = fullfile(Output_Dir, Output_File);
    
    if isfile(path)
        disp('FILE ALREADY EXISTS.')
    else
        % 使用 writematrix 写入数据
        writematrix(outputFeature, path);
    end
end


%}
