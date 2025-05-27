clear

%% load drugs date
load('TestPop\SR\female\SRfeaturesf.mat')
SR_feature = SRfeaturesf;

BaseFeName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'DCai', 'CTA', 'CTD50', 'CTD90', 'CTDtri', 'dCa'};

%% Amiodarone
%
load('TestPop\SRClassIII\female\Amiodarone\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr1 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Amiodarone\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dofetilide
%
load('TestPop\SRClassIII\female\Dofetilide\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr2 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Dofetilide\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dronedarone
%
load('TestPop\SRClassIII\female\Dronedarone\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr3 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Dronedarone\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Ibutilide
%
load('TestPop\SRClassIII\female\Ibutilide\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr4 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Ibutilide\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Sotalol
%
load('TestPop\SRClassIII\female\Sotalol\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr5 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Sotalol\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Vernakalant
%
load('TestPop\SRClassIII\female\Vernakalant\CABaseCells.mat') % population data7
[CABasefeture,Y_Arr] = Calibrate(SR_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr6 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\SRClassIII\female\Vernakalant\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

y = [rateArr1;rateArr2;rateArr3;rateArr4;rateArr5;rateArr6;30];

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
                idex = find(tx > 500) ;
                for i = 1:length(idex)
                    
                    if Vx(idex(i)) > (-30)
                         proArrhythmicity = 1;
                    end
                end
            end

end