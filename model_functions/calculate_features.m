function outputs = calculate_features(Vi,Caii,ti)

[t,V,Cai] = splitdata(Vi,Caii,ti); % saved the n beats, separate each beat
t = t - t(1); %为什么要减一下，计算微分           
Vderiv = diff(V)./diff(t) ;
[dVdtmax,dexmax] = max(Vderiv) ;
tinit = 100 ; %Time of maximum dV/dt, consider this beginning of action potential
vrest = min(V(1:dexmax));
[peakV,peakdex] = max(V) ;
tpeak = t(peakdex) ;
V20_exact = 0.8*(peakV - vrest) + vrest ;
V40_exact = 0.6*(peakV - vrest) + vrest ;
V50_exact = 0.5*(peakV - vrest) + vrest ;
V90_exact = 0.1*(peakV - vrest) + vrest ;
V20dex = find(t > tpeak & V < V20_exact) ;
V40dex = find(t > tpeak & V < V40_exact);
V50dex = find(t > tpeak & V < V50_exact) ;
V90dex = find(t > tpeak & V < V90_exact);
dCa = Cai(1);
DCai = max(Cai)-min(Cai);
[peakCa,peakdex] = max(Cai) ;
Capeak = peakCa;
tpeak_Ca = t(peakdex);
Ca90_exact = 0.1*(DCai) + Cai(end);
Ca90dex = find(t > tpeak_Ca & Cai < Ca90_exact);
Ca50_exact = 0.5*(DCai) + Cai(end) ;
Ca50dex = find(t > tpeak_Ca & Cai < Ca50_exact) ;

V20_time = t(V20dex(1));
APD20 = V20_time - tinit ;

APA = peakV - vrest ;
RMP= vrest;



if isempty(V40dex)
    APD40 = NaN;
else
    V40_time = t(V40dex(1));
    APD40 = V40_time - tinit;    
end

if isempty(V50dex)
    APD50 = NaN;
else
    V50_time = t(V50dex(1));
    APD50 = V50_time - tinit;    
end
    
if isempty(V90dex)    
    APD90 = NaN;
    TriAP = NaN;
else    
    V90_time = t(V90dex(1));
    APD90 = V90_time - tinit;             
    TriAP = V90_time-V50_time; %Triangulation
    
end


if isempty(Ca50dex)    
    CaD50 = NaN;
else    
    Ca50_time = t(Ca50dex(1));
    CaD50 = Ca50_time - tinit;    
end

if isempty(Ca90dex)
    CaD90 = NaN;
    TriCa = NaN;
else
    Ca90_time = t(Ca90dex(1));
    CaD90 = Ca90_time - tinit;
    TriCa = Ca90_time-Ca50_time;    
end
AP_metrics =[RMP dVdtmax APA APD20 APD40 APD50 APD90 TriAP];
CaT_metrics = [DCai Capeak CaD50 CaD90 TriCa dCa];

outputs = [AP_metrics CaT_metrics];
end

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
 function [times,volts,cais] = splitdata(Vi,Caii,Ti)   
            
            PCL = 1000;
            i=fix(length(Ti)./PCL);
            intervals = find(~mod(Ti,PCL));%查找ti取模为0的地方
            times = Ti(intervals(i):intervals(i+1),1);%分离已模拟的细胞的时间数据，用以后续检查数据的正确与否
            volts = Vi(intervals(i):intervals(i+1),1);%分离已模拟的细胞的电压数据
            cais = Caii(intervals(i):intervals(i+1),1);
           
 end