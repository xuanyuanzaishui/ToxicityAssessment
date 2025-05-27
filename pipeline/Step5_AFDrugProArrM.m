clear

%% load drugs date
load('TestPop\AF\male\AFfeaturesm.mat')
AF_feature = AFfeaturesm;

BaseFeName = {'RMP', 'dV/dtmax', 'APA', 'APD20', 'APD40', 'APD50',...
     'APD90','APDtri', 'DCai', 'CTA', 'CTD50', 'CTD90', 'CTDtri', 'dCa'};

%% Amiodarone
%
load('TestPop\AFClassIII\male\Amiodarone\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr1 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Amiodarone\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dofetilide
%
load('TestPop\AFClassIII\male\Dofetilide\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr2 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Dofetilide\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Dronedarone
%
load('TestPop\AFClassIII\male\Dronedarone\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr3 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Dronedarone\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Ibutilide
%
load('TestPop\AFClassIII\male\Ibutilide\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr4 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Ibutilide\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Sotalol
%
load('TestPop\AFClassIII\male\Sotalol\CABaseCells.mat') % population data
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr5 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Sotalol\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

%% Vernakalant
%
load('TestPop\AFClassIII\male\Vernakalant\CABaseCells.mat') % population data7
[CABasefeture,Y_Arr] = Calibrate(AF_feature,CABaseCells);

arrhy = 1;
for i = 1:length(Y_Arr)
    if Y_Arr(i)
        arrhy = arrhy + 1;
    end
end
rateArr6 = (arrhy-1)/length(Y_Arr);

matfile = fullfile('TestPop\AFClassIII\male\Vernakalant\', 'Y_Arr.mat');
save(matfile,'CABasefeture','Y_Arr')
%}

y = [rateArr1;rateArr2;rateArr3;rateArr4;rateArr5;rateArr6;length(Y_Arr)];

save('ym.mat','y')

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