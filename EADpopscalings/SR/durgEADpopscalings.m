clear

%% 1、查找EAD  2、保存scaling,用于python绘制分布图
%
%查找EAD male
load('TestPop/SR/male/SRpopICsm.mat','SRpopscalingsm');
load('TestPop/SR/female/SRpopICsf.mat','SRpopscalingsf');

settings.Folder = 'TestPop/SRClassIII/male';%SR
settings.Folder1 = 'TestPop/SRClassIII/female';%SR

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
        popscalingsm(arrhy,:) = SRpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 0;
popscalingsm(1:end,18) = 1;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = SRpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 0;
popscalingsf(1:end,18) = 1;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model2\conduction\'; % path
Output_File1 = 'EADpopscalingsmf2.csv';
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
        popscalingsm(arrhy,:) = SRpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 0;
popscalingsm(1:end,18) = 1;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = SRpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 0;
popscalingsf(1:end,18) = 1;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model2\conduction\'; % path
Output_File1 = 'EADpopscalingsmf3.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}


%% Ibutilide
%
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = SRpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 0;
popscalingsm(1:end,18) = 1;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = SRpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 0;
popscalingsf(1:end,18) = 1;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model2\conduction\'; % path
Output_File1 = 'EADpopscalingsmf4.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}

%% Sotalol
%
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = SRpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 0;
popscalingsm(1:end,18) = 1;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = SRpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 0;
popscalingsf(1:end,18) = 1;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model2\conduction\'; % path
Output_File1 = 'EADpopscalingsmf5.csv';
outputFile(Output_Dir,Output_File1,popscalingsmf)%特征
%}

%% Vernakalant
%
%male
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')
popscalingsm = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsm(arrhy,:) = SRpopscalingsm(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsm(1:end,15) = 0;
popscalingsm(1:end,16) = 1;
popscalingsm(1:end,17) = 0;
popscalingsm(1:end,18) = 1;

%female
yourFolder = fullfile(settings.Folder1,settings.SubFolder); 
File1 = fullfile(yourFolder,'Y_Arr.mat');
load(File1,'Y_Arr')

popscalingsf = [];
arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        popscalingsf(arrhy,:) = SRpopscalingsf(i,:);
        arrhy = arrhy + 1;
    end
end
popscalingsf(1:end,15) = 1;
popscalingsf(1:end,16) = 0;
popscalingsf(1:end,17) = 0;
popscalingsf(1:end,18) = 1;

popscalingsmf = [popscalingsm;popscalingsf];
Output_Dir = 'E:\matlab\model2\conduction\'; % path
Output_File1 = 'EADpopscalingsmf6.csv';
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
