function APD90 = find_APD90(t,V)
%--------------------------------------------------------------------------
                        %% -- find_APD90.m -- %%
% Description: computes action potential duration at 90% repolarization 

% Inputs:
% --> t - [double array] time array   
% --> V - [double array] voltage array 

% Outputs: 
% --> APD90 - [double] computes action potential duration at 90%
% repolarization
% -------------------------------------------------------------------------
%%
t = t - t(1); %为什么要减一下，计算微分
Vderiv = diff(V)./diff(t) ;
[~,dexmax] = max(Vderiv) ;%返回最大值的索引
tinit = t(dexmax(1));  %Time of maximum dV/dt, consider this beginning of action potential
vrest = min(V(1:dexmax));%返回静息电位的值
[peakV,peakdex] = max(V) ;%返回电压最大值及其索引
tpeak = t(peakdex) ;%返回电压最大值对应的时间
V90_exact = (0.1*(peakV - vrest) + vrest) ;

V90dex = find(t > tpeak & V < V90_exact);

if ~isempty(V90dex)
    V90_time = t(V90dex(1));
    APD90 = V90_time - 100;
else
    APD90 = NaN;
end
