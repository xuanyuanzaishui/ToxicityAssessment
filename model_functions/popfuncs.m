classdef popfuncs
    
    methods(Static)
%--------------------------------------------------------------------------
                      %% -- create_scale_vector.m -- %%
% Description: creates population variability scaling matrix.  

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> variations - [double array] number of members of the population 

% Outputs: 
% --> scaling_matrix - [double array] population variability scaling matrix
% -------------------------------------------------------------------------
        function scaling_matrix = create_scale_vector(settings,variations)                     
            [~,c] = model_parameters('endo');
            nG = length(fieldnames(c.G));
            
            scaling_matrix = exp(settings.sigmaG*randn(nG,variations))' ; % same parameter variation for each pop
           
        end 
               
%--------------------------------------------------------------------------
                 %% -- clean_population.m -- %%
% Description: Determines if there are arrhythmias in the population,
% removes them, and reruns new cells 

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Xin - [struct] population time and state variables  

% Outputs: 
% --> Xout - [struct] cleaned population 
% -------------------------------------------------------------------------
        function Xout = clean_population(Xin,settings)
            flag = 1; % run the while loop
            settings_rerun = settings;
            X2 = Xin; 
            total_variations = settings.variations;

            while flag     %生成一个种群设置细胞的数目，该种群的细胞数据均是符合要求的 
                if settings.remove_arrhythmias
                    removes_1 = popfuncs.find_arrhythmias(settings_rerun,X2); 
                                        % remove cells with arrhythmias before trigger
                    X2(removes_1) = [];%清除发生失常的细胞数据
                end
                n_X2 = length(X2);  % #of cells kept
                
                if length(X2) == total_variations || settings.reruncells == 0 
                                                   %如果数目不够，决定是否要重新模拟
                    flag = 0;
                else
                    settings_rerun.variations = total_variations - n_X2; % number of cells to rerun
                    settings_rerun.scalings = popfuncs.create_scale_vector(settings_rerun,settings_rerun.variations);
                    disp(['Need to Rerun: ' num2str(settings_rerun.variations) ' cells.'])
                    
                    pert = settings_blockcurrents;
                    X = runSim(settings_rerun,pert); % run population simulation%再次生成模拟细胞
                    X2 = [X2,X];
                end
            end
            Xout = X2;
        end        
%--------------------------------------------------------------------------
  
                 %% -- find_arrhythmias.m -- %%失常数据
% Description: Determines if there are arrhythmias in a population

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Xin - [struct] population time and state variables  

% Outputs: 
% --> remove_AP - [logical] logical array of which cells have arrhythmias  
% -------------------------------------------------------------------------
  function remove_AP = find_arrhythmias(settings,Xin)
           remove_AP = false(length(Xin),1);
        for ii = 1:length(Xin) %for each member of the population 
            Xt = Xin(ii).times;
            XStates = Xin(ii).states;
            [times,volts] = popfuncs.splitdata(Xt,XStates,settings); 
                             % saved the last 10 beats, separate each beat
           % check each beat for arrhythmic activity 
           for i = 1:settings.numbertokeep %for each beat 检查保存的节拍数据
               t = times{i}; 
               t = t - t(1);
               V = volts{i}(:,1);
               APDs = find_APD90(t,V);%加载APD的计算函数文件

               %修改位置
              if isnan(APDs) || V(1) > -63 || max(V) < 0 || V(end)> -63 %check if cell failed to repolarize
                 remove_AP(ii) = true;                      %检查细胞是否重新极化失败，极化失败则移除
                                                               
              else %check AP for early afterdepolarizations否则就检查动作电位的早期后去极化是否正确
                   [~,index] = find(t==settings.stim_delay); %找到有12.5刺激开始的位置
                   time_roi = t(index:end); %从延迟开始到模拟结束的时间序列
                   Vm_roi = V(index:end);   %从延迟开始到模拟结束的电压序列
                   dVm_roi = (Vm_roi(2:end)-Vm_roi(1:end-1))./(time_roi(2:end)-time_roi(1:end-1));%电压对时间的微分
                   dVm_dVm = dVm_roi(1:end-1).*dVm_roi(2:end);
                   [~, idx_dVm] = max(dVm_roi);%找到最大微分的位置
                   dVm_0 = dVm_dVm(idx_dVm:end) < -1.0*10^-6; % -1.0*10^-6存在问题0或1
                   dVm_t = time_roi(3:end);
                   tVm_0 = dVm_t(idx_dVm:end);
                   ts = tVm_0(dVm_0);
                   time_between_depols = diff(ts);
                        
                 if any(time_between_depols > 130) %如果  的值大于130则移除
                      remove_AP(ii) = true;  
                 end  
              end
              if APDs < 251 && APDs > 385%check for alternans SR 
                   remove_AP(ii) = true;
              end
           end
        end
  end
%--------------------------------------------------------------------------
                 %% -- calibrate_experimental.m -- 
% Description:

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Ti - [double array] time matrix 
% --> V - [double array] voltage matrix 

% Outputs: 
% --> times - [cell] time vector for each beat 
% --> volts - [cell] voltage vector for each beat  
%--------------------------------------------------------------------------
function calibrate = calibrate_experimental(settings,Xin,tflage)
                X2=[];
                Y2=[];
                Z2=[];
       for ii = 1:length(Xin) %for each member of the population
           Xt = Xin(ii).times;
           XStates = Xin(ii).states;
           [times,volts,cais] = popfuncs.splitdata(Xt,XStates,settings); 
                                    % saved the last 10 beats, separate each beat
           t = times{end};
           V = volts{end};
           Cai = cais{end} * 1000000;
           %%%预计修改位置
           outputs = calculate_features(V,Cai,t);%加载特征计算函数文件
           RMP = outputs(1,1);
           dVdt = outputs(1,2);
           APA = outputs(1,3);
           APD20 = outputs(1,4);
           APD50 = outputs(1,6);
           APD90 = outputs(1,7);
           DCai = outputs(1,9);
           switch tflage
                case 1
                    I1 = RMP > (-85) && RMP < (-65);
                    I2 = APA > 80 && APA < 130;
                    I3 = APD20 > 1 && APD20 < 75;%定义误差范围，判断模拟数据是否正常
                    I4 = APD50 > 30 && APD50 < 180;
                    I5 = APD90 > 190 && APD90 < 330;
                    I6 = dVdt > 40 && dVdt < 420;
                   
                    I = [I1,I2,I3,I4,I5,I6];
                    if all(I == 1)
                      X1=Xin(ii);
                      X2 = [X2,X1];
                    end

               case 2
                    J1 = RMP > (-85.2) && RMP < (-66.1);
                    J2 = APA > 94.6 && APA < 137.5;
                    J3 = APD20 > 0.7 && APD20 < 48;%定义误差范围，判断模拟数据是否正常
                    J4 = APD50 > 30 && APD50 < 180;
                    J5 = APD90 > 140 && APD90 < 330;
                    J6 = dVdt > 40 && dVdt < 420;
                    J7 = DCai > 100 ;

                    J = [J1,J2,J3,J4,J5,J6,J7];
                    if all(J == 1)
                       Y1=Xin(ii);
                       Y2 = [Y2,Y1];
                    end

               case 3
                    K1 = RMP > (-90) && RMP < (-66);
                    K2 = APA > 102.6 && APA < 130.6;
                    K3 = APD20 > 12 && APD20 < 48;%定义误差范围，判断模拟数据是否正常
                    K4 = APD50 > 35.2 && APD50 < 109.2;
                    K5 = APD90 > 138 && APD90 < 262;

                    K = [K1,K2,K3,K4,K5];
                    if all(K == 1)
                      Z1=Xin(ii);
                      Z2 = [Z2,Z1];
                    end
           end
       end
       switch tflage
            case 1
                  calibrate=X2;
            case 2
                  calibrate=Y2;
            case 3
                  calibrate=Z2; 
       end
end
%--------------------------------------------------------------------------
                      %% -- splitdata.m -- %%
% Description: When multiple beats of a single cell are simulated, this
% function separates each beat into its own cell array. This is used
% mainly when settings.numbertokeep is greater than 1. 

% Inputs:
% --> settings - [struct] simulation protocol settings    
% --> Ti - [double array] time matrix 
% --> V - [double array] voltage matrix 

% Outputs: 
% --> times - [cell] time vector for each beat 
% --> volts - [cell] voltage vector for each beat 
% -------------------------------------------------------------------------
        function [times,volts,cais] = splitdata(Ti,States,settings)   
            numberofAPs = settings.numbertokeep;
            PCL = settings.PCL;
            intervals = find(~mod(Ti,PCL));%查找ti取模为0的地方
            times = {};
            volts ={};
            cais = {};
            for i = 1:numberofAPs
                times{end+1} = Ti(intervals(i):intervals(i+1));%分离已模拟的细胞的时间数据，用以后续检查数据的正确与否
                volts{end+1} = States(intervals(i):intervals(i+1),39);%分离已模拟的细胞的电压数据
                cais{end+1} = States(intervals(i):intervals(i+1),38);%分离已模拟的细胞的钙离子数据
            end
        end
    end 
end
