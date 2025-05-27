function datatable = runSim(settings,pert)
%% 1--- Load Model Parameters 
[p,c] = model_parameters(settings.celltype);%加载函数文件

%% 2--- Load Model Initial Conditions of State Variables 加载参数初始值文件
if ~isfield(settings,'ICs')
    y0 = ICs(settings.steady_state,settings.PCL);
    y0s = repmat(y0,settings.variations,1);   %repmat复制和平铺矩阵
else 
    y0s = settings.ICs;
end

%% 3--- Define Simulation Protocol 
stim_starts = settings.stim_delay + settings.PCL*(0:settings.nBeats-1) ;%开始定义
stim_ends = stim_starts + settings.stim_dur ;%结束定义

% Create intervals for each beat 刺激间隔1000ms
simints = 3*settings.nBeats ; %settings.nBeats施加的刺激数量
for i=1:settings.nBeats
    intervals(3*i-2,:) = [settings.PCL*(i-1),stim_starts(i)] ; %beginning 0-100
    intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ; %stimulus 极短时间内施加刺激100-105
    intervals(3*i,:) = [stim_ends(i),settings.PCL*i] ; % stimulus ends 105-1100
end
tend = settings.nBeats*settings.PCL ;              % end of simulation, ms
intervals(end,:) = [stim_ends(end),tend] ;

% Determine when to apply stim_amp or 0 amp 什么时候施加刺激
Istim = zeros(simints,1) ;
stimindices = 3*(1:settings.nBeats) - 1 ; % apply stimulus on second part of intervals
Istim(stimindices) = settings.stim_amp ; 

%% 4--- Population Variables 
F_G = fieldnames(c.G);
S = settings.scalings;

%% 5--- Define Perturbation Protocol 
[p,c]= perturbations(c,p,pert);  %关键在于c中数值的改变，p不变
baselineparameters = c;

%% 6--- Run Loop 
for ii=1:settings.variations    
    scaling = S(ii,:);          %加载正态分布数值
    c = scaling_factors(scaling,baselineparameters,F_G); %主要是通道在变
    statevar_i = y0s(ii,:);     %加载参数初始值
    options = odeset('RelTol',1e-5,'MaxStep',1);
    
    % stimulate cell
    t = 0 ;
    statevars = statevar_i ;
    for i=1:simints            %一个细胞节拍数目
        [post,posstatevars] = ode15s(@dydt_Ohara,intervals(i,:),statevar_i,options,Istim(i),p,c);
         
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)];
        statevar_i = posstatevars(end,:);
    end 
    % Only save the number of beats specified in numbertokeep
    start = find( t == intervals(simints-3*settings.numbertokeep+1,1));
    t_final = t(start:end);    
    statevars_final = statevars(start:end,:);

    datatable(ii).times =  t_final - min(t_final) ;
     datatable(ii).states = statevars_final;
    datatable(ii).scalings = scaling;
    %
    currents=[];
    for i=1:size(t_final)
        currents=[currents;dydt_Ohara(t_final(i),statevars_final(i,:),20,p,c,0)];
    end
    datatable(ii).currents = currents;
    %}                         
end

%% Nested Functions 
function [p,c]= perturbations(c,p,pert)
    c.G.GNa  = c.G.GNa * pert.GNa;
    c.G.GNaL = c.G.GNaL * pert.GNaL;
    c.G.Gto  = c.G.Gto * pert.Gto;
    c.G.GKr = c.G.GKr * pert.GKr;
    c.G.GKs = c.G.GKs * pert.GKs;
    c.G.GKur = c.G.GKur * pert.GKur;
    c.G.GK1  = c.G.GK1 * pert.GK1;
    c.G.GCaL = c.G.GCaL * pert.GCaL;
    c.G.Pnak = c.G.Pnak * pert.NaK;
    c.G.Gncx = c.G.Gncx * pert.GNCX;
    c.G.GpCa = c.G.GpCa * pert.GpCa;
    c.G.Rel   = c.G.Rel * pert.Rel ;
    c.G.SERCA = c.G.SERCA * pert.SERCA ;
    c.G.Leak  = c.G.Leak * pert.Leak;


function c = scaling_factors(scaling,baselineparameters,n_G)
    scalings.G = scaling(1:length(n_G));
    for iF = 1:length(n_G)
        aF = n_G{iF};
        c.G.(aF) = baselineparameters.G.(aF) * scalings.G(iF);
    end 