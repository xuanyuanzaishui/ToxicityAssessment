clear

%% load drugs date
load('TestPop\AF\female\AFfeaturesf.mat')
AF_feature = AFfeaturesf;

BaseFeName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'DCai', 'CTA', 'CTD50', 'CTD90', 'CTDtri', 'dCa'};


settings.Folder = 'TestPop/AFClassIII/female/';%SR
%% Amiodarone
%
settings.SubFolder = 'Amiodarone';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr1 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dofetilide
%
settings.SubFolder = 'Dofetilide';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr2 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dronedarone
%
settings.SubFolder = 'Dronedarone';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr3 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Ibutilide
%
settings.SubFolder = 'Ibutilide';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr4 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Sotalol
%
settings.SubFolder = 'Sotalol';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr5 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Vernakalant
%
settings.SubFolder = 'Vernakalant';
yourFolder = fullfile(settings.Folder,settings.SubFolder); 
File = fullfile(yourFolder,'CABaseCells.mat'); 
load(File,'CABaseCells');
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr6 = (arrhy-1)/length(Y_Arr);

matfile = fullfile(yourFolder, 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

y = [rateArr1;rateArr2;rateArr3;rateArr4;rateArr5;rateArr6;length(Y_Arr)];

save('yf.mat','y')

%% 函数

%促心律失常
function [CABaseFea,Y_Arr] = Calibrate(AF_feature,drugBaseCells)
            CABaseFea =[];
            Y_Arr = zeros(length(AF_feature),1);
            
            for i = 1:length(drugBaseCells)
                ti = drugBaseCells(i).times;
                Vi = drugBaseCells(i).V;
                Caii = drugBaseCells(i).Cai;
                [t,V,Cai] = splitdata(Vi,Caii,ti);

                CABaseFea(i,:) = calculate_features(V,Cai,t);      
                
                arrhythmicity = proArr(t,V,CABaseFea(i,:));

                if arrhythmicity
                    Y_Arr(i,1) = 1;
                end
        
            end
end

%数据拆分
function [times,volts,cais] = splitdata(Vi,Caii,Ti)   
            PCL = 1000;
            i=fix(length(Ti)./PCL);
            intervals = find(~mod(Ti,PCL));%查找ti取模为0的地方
            times = Ti(intervals(i):intervals(i+1),1);%分离已模拟的细胞的时间数据，用以后续检查数据的正确与否
            volts = Vi(intervals(i):intervals(i+1),1);%分离已模拟的细胞的电压数据
            cais = Caii(intervals(i):intervals(i+1),1);
           
end

%
function proArrhythmicity = proArr(tx,Vx,CABaseFeature)
           proArrhythmicity = 0;
            if anynan(CABaseFeature)
                proArrhythmicity = 1;
            else
                tx = tx - tx(1);       
                Vderiv = diff(Vx)./diff(tx);
                [~,peakdex] = max(Vx) ;
                idex = find(tx > tx(peakdex) & Vx < -14) ;
                for i = 520 : 1000
                    if Vderiv(i)>0
                        proArrhythmicity = 1;
                    end
                end
            
            end

end