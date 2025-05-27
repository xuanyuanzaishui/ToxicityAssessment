function out = dydt_Grandi(t,statevar,Id,~,c,flag)
if ~exist('flag','var') || flag 
    flag = 1;
end 

yxt=statevar';
ydot = zeros(size(yxt));

%% Model Parameters
AF=0;
% ISO
ISO=0;
% Right ATRIUM
RA=1;
% Constants
R = 8314.0;      %[J/kmol*K]  
Frdy = 96485.0;    %[C/mol]      
Temp = 310.0;     % [K]       
FoRT = Frdy/R/Temp;
Cmem = 1.1e-10;   % [F] membrane capacitance 1.3810e-10;
Qpow = (Temp-310)/10;

% Cell geometry几何
cellLength = 100;     % cell length [um]113;%100 
cellRadius = 10.25;   % cell radius [um]12;%10.25
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]

Vcell = pi*cellRadius^2*cellLength*1e-15;  % [L]
Vmyo = 0.65*Vcell;
Vsr = 0.035*Vcell; 
Vsl = 0.02*Vcell; 
Vjunc = 1*0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;% [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2] 

J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments间中的分数电流
Fjunc = 0.11;   
Fsl = 1-Fjunc;
Fjunc_CaL = 0.9;
Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations固定离子浓度    
Cli = 15;   % Intracellular Cl  [mM]xibaonei
Clo = 150;  % Extracellular Cl  [mM]xibaowai
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/yxt(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/yxt(33));       % [mV]
ek = (1/FoRT)*log(Ko/yxt(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/yxt(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/yxt(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

% Na currents/transport parameters
GNa=c.G.GNa;  % [mS/uF] [c.G.GNa ]
GNaB = 0.597e-3;      % [mS/uF]  [c.G.PNab ]
IbarNaK = c.G.Pnak;     % [uA/uF]
KmNaip = 11*(1-0.25*ISO);   % [mM]
KmKo =1.5;         % [mM]

% K current parameters
pNaK = 0.01833;   %  
gkp = 0.002;    %[c.G.GKb 37]

% Cl current parameters
GClCa =0.0548;    % [mS/uF]
GClB = 9e-3;       % [mS/uF]
KdClCa = 100e-3;    % [mM]

% I_Ca parameters
pNa = (1+0.5*ISO)*0.75e-8*c.G.GCaL;       % [cm/sec]
pCa = (1+0.5*ISO)*2.7e-4*c.G.GCaL;       % [cm/sec]
pK = (1+0.5*ISO)*1.35e-7*c.G.GCaL;        % [cm/sec]
Q10CaL = 1.8;  

% Ca transport parameters
IbarNCX = c.G.Gncx; %[c.G.Gncx 36]     % [uA/uF]5.5 before - 9 in rabbit
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]  
nu = 0.35;          % [none]
Kdact =0.384e-3;   % [mM] 0.256 rabbit
Q10NCX = 1.57;      % [none]
IbarSLCaP = c.G.GpCa; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa =0.5e-3;     % [mM] [c.G.GpCa 38]
GCaB = 6.0643e-4;    % [uA/uF] [c.G.PCab 43]
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = 5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = (2.5-1.25*ISO)*0.246e-3;          % [mM] default
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 25;                 % [1/ms]      
koCa = 10+20*AF+10*ISO*(1-AF);               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = (1+0.5*ISO)*19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

%% Membrane Currents
mss = 1 / ((1 + exp( -(56.86 + yxt(39)) / 9.03 ))^2);
taum = 0.1292 * exp(-((yxt(39)+45.79)/15.54)^2) + 0.06487 * exp(-((yxt(39)-4.823)/51.12)^2);                 
 
ah = (yxt(39) >= -40) * (0)... 
   + (yxt(39) < -40) * (0.057 * exp( -(yxt(39) + 80) / 6.8 )); 
bh = (yxt(39) >= -40) * (0.77 / (0.13*(1 + exp( -(yxt(39) + 10.66) / 11.1 )))) ...
   + (yxt(39) < -40) * ((2.7 * exp( 0.079 * yxt(39)) + 3.1*10^5 * exp(0.3485 * yxt(39)))); 
tauh = 1 / (ah + bh); 
hss = 1 / ((1 + exp( (yxt(39) + 71.55)/7.43 ))^2);
 
aj = (yxt(39) >= -40) * (0) ...
    +(yxt(39) < -40) * (((-2.5428 * 10^4*exp(0.2444*yxt(39)) - 6.948*10^-6 * exp(-0.04391*yxt(39))) * (yxt(39) + 37.78)) / ...
                     (1 + exp( 0.311 * (yxt(39) + 79.23) )));
bj = (yxt(39) >= -40) * ((0.6 * exp( 0.057 * yxt(39))) / (1 + exp( -0.1 * (yxt(39) + 32) ))) ...
   + (yxt(39) < -40) * ((0.02424 * exp( -0.01052 * yxt(39) )) / (1 + exp( -0.1378 * (yxt(39) + 40.14) ))); 
tauj = 1 / (aj + bj);
jss = 1 / ((1 + exp( (yxt(39) + 71.55)/7.43 ))^2);         
 
ydot(1) = (mss - yxt(1)) / taum;
ydot(2) = (hss - yxt(2)) / tauh;
ydot(3) = (jss - yxt(3)) / tauj;
    
I_Na_junc = Fjunc*GNa*yxt(1)^3*yxt(2)*yxt(3)*(yxt(39)-ena_junc);
I_Na_sl = Fsl*GNa*yxt(1)^3*yxt(2)*yxt(3)*(yxt(39)-ena_sl);
I_Na = I_Na_junc+I_Na_sl;

% Late I_Na
GNaL=c.G.GNaL; %[c.G.GNaL 31]
aml = 0.32*(yxt(39)+47.13)/(1-exp(-0.1*(yxt(39)+47.13)));
bml = 0.08*exp(-yxt(39)/11);
hlinf = 1/(1+exp((yxt(39)+91)/6.1));
tauhl=600;
ydot(60) = aml*(1-yxt(60))-bml*yxt(60);
ydot(61) = (hlinf-yxt(61))/tauhl;

I_NaL_junc = Fjunc*GNaL*yxt(60)^3*yxt(61)*(yxt(39)-ena_junc);
I_NaL_sl = Fsl*GNaL*yxt(60)^3*yxt(61)*(yxt(39)-ena_sl);
I_NaL = I_NaL_junc + I_NaL_sl;
if t<9050
    ydot(62)=0;
else
    ydot(62)=I_NaL;
end
% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(yxt(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(yxt(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*yxt(39)*FoRT)+0.0365*sigma*exp(-yxt(39)*FoRT));
I_nak_junc = 1*Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/yxt(32))^4) /(Ko+KmKo);
I_nak_sl = 1*Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/yxt(33))^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
gkr=(c.G.GKr)*sqrt(Ko/5.4); %[c.G.GKr_ 33]
xrss = 1/(1+exp(-(yxt(39)+10)/5));
tauxr = 550/(1+exp((-22-yxt(39))/9))*6/(1+exp((yxt(39)-(-11))/9))+230/(1+exp((yxt(39)-(-40))/20));
ydot(12) = (xrss-yxt(12))/tauxr;
rkr = 1/(1+exp((yxt(39)+74)/24));
I_kr = gkr*yxt(12)*rkr*(yxt(39)-ek);


%% I_ks: Slowly Activating K Current
markov_iks=0;
% pcaks_junc = -log10(y(36))+3.0; 
% pcaks_sl = -log10(y(37))+3.0;  
% gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
% gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     

eks = (1/FoRT)*log((Ko+pNaK*Nao)/(yxt(35)+pNaK*yxt(34)));

if markov_iks==0
gks_junc=(c.G.GKs)*1; %[c.G.GKs_ 34]
gks_sl=(c.G.GKs)*1;  %[c.G.GKs_ 34]
xsss = 1 / (1+exp(-(yxt(39)+40*ISO + 3.8)/14.25)); % fitting Fra
tauxs=990.1/(1+exp(-(yxt(39)+40*ISO+2.436)/14.12));
ydot(13) = (xsss-yxt(13))/tauxs;
I_ks_junc = Fjunc*gks_junc*yxt(13)^2*(yxt(39)-eks);
I_ks_sl = Fsl*gks_sl*yxt(13)^2*(yxt(39)-eks);                                                                                                                                   
I_ks = I_ks_junc+I_ks_sl;
else
    gks_junc=1*0.0065;
    gks_sl=1*0.0065; %FRA
    alpha=3.98e-4*exp(3.61e-1*yxt(39)*FoRT);
    beta=5.74e-5*exp(-9.23e-2*yxt(39)*FoRT);
    gamma=3.41e-3*exp(8.68e-1*yxt(39)*FoRT);
    delta=1.2e-3*exp(-3.3e-1*yxt(39)*FoRT);
    teta=6.47e-3;
    eta=1.25e-2*exp(-4.81e-1*yxt(39)*FoRT);
    psi=6.33e-3*exp(1.27*yxt(39)*FoRT);
    omega=4.91e-3*exp(-6.79e-1*yxt(39)*FoRT);
    
    ydot(42)=-4*alpha*yxt(42)+beta*yxt(43);
    ydot(43)=4*alpha*yxt(42)-(beta+gamma+3*alpha)*yxt(43)+2*beta*yxt(44);
    ydot(44)=3*alpha*yxt(43)-(2*beta+2*gamma+2*alpha)*yxt(44)+3*beta*yxt(45);
    ydot(45)=2*alpha*yxt(44)-(3*beta+3*gamma+alpha)*yxt(45)+4*beta*yxt(46);
    ydot(46)=1*alpha*yxt(44)-(4*beta+4*gamma)*yxt(46)+delta*yxt(50);    
    ydot(47)=gamma*yxt(43)-(delta+3*alpha)*yxt(47)+beta*yxt(48);   
    ydot(48)=2*gamma*yxt(44)+3*alpha*yxt(47)-(delta+beta+2*alpha+gamma)*yxt(48)+2*beta*yxt(49)+2*delta*yxt(51);
    ydot(49)=3*gamma*yxt(45)+2*alpha*yxt(48)-(delta+2*beta+1*alpha+2*gamma)*yxt(49)+3*beta*yxt(50)+2*delta*yxt(52);
    ydot(50)=4*gamma*yxt(46)+1*alpha*yxt(49)-(delta+3*beta+0*alpha+3*gamma)*yxt(50)+2*delta*yxt(53);
    ydot(51)=1*gamma*yxt(48)-(2*delta+2*alpha)*yxt(51)+beta*yxt(52);  
    ydot(52)=2*gamma*yxt(49)+2*alpha*yxt(51)-(2*delta+beta+1*alpha+gamma)*yxt(52)+2*beta*yxt(53)+3*delta*yxt(54);
    ydot(53)=3*gamma*yxt(50)+1*alpha*yxt(52)-(2*delta+2*beta+2*gamma)*yxt(53)+3*delta*yxt(55);
    ydot(54)=1*gamma*yxt(52)-(3*delta+1*alpha)*yxt(54)+beta*yxt(55);  
    ydot(55)=2*gamma*yxt(53)+1*alpha*yxt(54)-(3*delta+1*beta+1*gamma)*yxt(55)+4*delta*yxt(56);
    ydot(56)=1*gamma*yxt(55)-(4*delta+teta)*yxt(56)+eta*yxt(57);
    O2=1-(yxt(42)+yxt(43)+yxt(44)+yxt(45)+yxt(46)+yxt(47)+yxt(49)+yxt(48)+yxt(50)+yxt(51)+yxt(52)+yxt(53)+yxt(54)+yxt(55)+yxt(56)+yxt(57));
    ydot(57)=1*teta*yxt(56)-(eta+psi)*yxt(57)+omega*O2;
    I_ks_junc = Fjunc*gks_junc*(yxt(57)+O2)*(yxt(39)-eks);
    I_ks_sl = Fsl*gks_sl*(yxt(57)+O2)*(yxt(39)-eks);                                                                                                                                   
    I_ks = I_ks_junc+I_ks_sl;
end
%I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-yxt(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(yxt(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(yxt(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes

GtoFast=c.G.Gto; %nS/pF maleckar; %human atrium[c.G.Gto 32]

%11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
%atrium
%equations for activation; 
xtoss = ( (1)./ ( 1 + exp( -(yxt(39)+1.0)/11.0 ) ) );
tauxtof = 3.5*exp(-((yxt(39)/30.0)^2.0))+1.5;

%equations for inactivation;
ytoss = ( (1.0)./ ( 1 + exp( (yxt(39)+40.5)/11.5) ) ) ;
tauytof =25.635*exp(-(((yxt(39)+52.45)/15.8827)^2.0))+24.14;%14.14

ydot(10) = (xtoss-yxt(10))/tauxtof;
ydot(11) = (ytoss-yxt(11))/tauytof;
I_tof = 1.0*GtoFast*yxt(10)*yxt(11)*(yxt(39)-ek);
I_tos = 0;
I_to = I_tos + I_tof;

%% I_kur: Ultra rapid delayed rectifier Outward K Current
%Equation for IKur; from Maleckar et al. 2009 - EG
%atrium
%equations for activation;
Gkur = c.G.GKur*(1+0.2*RA); %nS/pF maleckar[c.G.GKur 39] 0.045
xkurss = ( (1)./ ( 1 + exp( (yxt(39)+6)/-8.6 ) ) );
tauxkur = 9/(1+exp((yxt(39)+5)/12.0))+0.5;

%equations for inactivation;
ykurss = ( (1)./ ( 1 + exp( (yxt(39)+7.5)/10 ) ) );
tauykur = 590/(1+exp((yxt(39)+60)/10.0))+3050;

ydot(58) = (xkurss-yxt(58))/tauxkur;
ydot(59) = (ykurss-yxt(59))/tauykur;
I_kur = 1*Gkur*yxt(58)*yxt(59)*(yxt(39)-ek);

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(yxt(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(yxt(39)+5.476-ek)) + exp(0.06175*(yxt(39)-ek-594.31))) /(1 + exp(-0.5143*(yxt(39)-ek+4.753)));
kiss = aki/(aki+bki);

%I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(y(39)-ek);
%SVP 11/11/09
%multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
%resting potential将IK1乘以0.15，将其扩展到单细胞分离的心房细胞静息电位
I_ki =c.G.GK1*sqrt(Ko/5.4)*kiss*(yxt(39)-ek);

% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/yxt(36))*(yxt(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/yxt(37))*(yxt(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(yxt(39)-ecl);

GClCFTR=0;%4.9e-3*ISO;     % [mS/uF]
I_ClCFTR = GClCFTR*(yxt(39)-ecl);

%% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(yxt(39)+3*ISO+9)/6)); %in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
taud = 1*dss*(1-exp(-(yxt(39)+3*ISO+9)/6))/(0.035*(yxt(39)+3*ISO+9)); 
fss = 1/(1+exp((yxt(39)+3*ISO+30)/7))+0.2/(1+exp((50-yxt(39)-3*ISO)/20)); % in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
tauf = 1/(0.0197*exp( -(0.0337*(yxt(39)+3*ISO+25))^2 )+0.02);
ydot(4) = (dss-yxt(4))/taud;
ydot(5) = (fss-yxt(5))/tauf;
ydot(6) = 1.7*yxt(36)*(1-yxt(6))-1*11.9e-3*yxt(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*yxt(37)*(1-yxt(7))-1*11.9e-3*yxt(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/yxt(37)));
fcaCaj= 0.1/(1+(0.01/yxt(36)));
fcaCaMSL=0;
fcaCaj= 0;
ibarca_j = pCa*4*(yxt(39)*Frdy*FoRT) * (0.341*yxt(36)*exp(2*yxt(39)*FoRT)-0.341*Cao) /(exp(2*yxt(39)*FoRT)-1);
ibarca_sl = pCa*4*(yxt(39)*Frdy*FoRT) * (0.341*yxt(37)*exp(2*yxt(39)*FoRT)-0.341*Cao) /(exp(2*yxt(39)*FoRT)-1);
ibark = pK*(yxt(39)*Frdy*FoRT)*(0.75*yxt(35)*exp(yxt(39)*FoRT)-0.75*Ko) /(exp(yxt(39)*FoRT)-1);
ibarna_j = pNa*(yxt(39)*Frdy*FoRT) *(0.75*yxt(32)*exp(yxt(39)*FoRT)-0.75*Nao)  /(exp(yxt(39)*FoRT)-1);
ibarna_sl = pNa*(yxt(39)*Frdy*FoRT) *(0.75*yxt(33)*exp(yxt(39)*FoRT)-0.75*Nao)  /(exp(yxt(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*yxt(4)*yxt(5)*((1-yxt(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_Ca_sl = (Fsl_CaL*ibarca_sl*yxt(4)*yxt(5)*((1-yxt(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*yxt(4)*yxt(5)*(Fjunc_CaL*(fcaCaj+(1-yxt(6)))+Fsl_CaL*(fcaCaMSL+(1-yxt(7))))*Q10CaL^Qpow)*0.45;
I_CaNa_junc = (Fjunc_CaL*ibarna_j*yxt(4)*yxt(5)*((1-yxt(6))+fcaCaj)*Q10CaL^Qpow)*0.45;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*yxt(4)*yxt(5)*((1-yxt(7))+fcaCaMSL)*Q10CaL^Qpow)*.45;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

% I_ncx: Na/Ca Exchanger flux
Ka_junc = 1/(1+(Kdact/yxt(36))^2);
Ka_sl = 1/(1+(Kdact/yxt(37))^2);
s1_junc = exp(nu*yxt(39)*FoRT)*yxt(32)^3*Cao;
s1_sl = exp(nu*yxt(39)*FoRT)*yxt(33)^3*Cao;
s2_junc = exp((nu-1)*yxt(39)*FoRT)*Nao^3*yxt(36);
s3_junc = KmCai*Nao^3*(1+(yxt(32)/KmNai)^3) + KmNao^3*yxt(36)*(1+yxt(36)/KmCai)+KmCao*yxt(32)^3+yxt(32)^3*Cao+Nao^3*yxt(36);
s2_sl = exp((nu-1)*yxt(39)*FoRT)*Nao^3*yxt(37);
s3_sl = KmCai*Nao^3*(1+(yxt(33)/KmNai)^3) + KmNao^3*yxt(37)*(1+yxt(37)/KmCai)+KmCao*yxt(33)^3+yxt(33)^3*Cao+Nao^3*yxt(37);

I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*yxt(39)*FoRT));
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*yxt(39)*FoRT));
I_ncx = I_ncx_junc+I_ncx_sl;

% I_pca: Sarcolemmal Ca Pump Current
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*yxt(36)^1.6/(KmPCa^1.6+yxt(36)^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*yxt(37)^1.6/(KmPCa^1.6+yxt(37)^1.6);
I_pca = I_pca_junc+I_pca_sl;

% I_cabk: Ca Background Current
I_cabk_junc = Fjunc*GCaB*(yxt(39)-eca_junc);
I_cabk_sl = Fsl*GCaB*(yxt(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/yxt(31))^2.5);
koSRCa = (1)*koCa/kCaSR;%
kiSRCa = kiCa*kCaSR;
RI = 1-yxt(14)-yxt(15)-yxt(16);
ydot(14) = (kim*RI-kiSRCa*yxt(36)*yxt(14))-(koSRCa*yxt(36)^2*yxt(14)-kom*yxt(15));   % R
ydot(15) = (koSRCa*yxt(36)^2*yxt(14)-kom*yxt(15))-(kiSRCa*yxt(36)*yxt(15)-kim*yxt(16));% O
ydot(16) = (kiSRCa*yxt(36)*yxt(15)-kim*yxt(16))-(kom*yxt(16)-koSRCa*yxt(36)^2*RI);   % I
J_SRCarel = ks*yxt(15)*(yxt(31)-yxt(36))*c.G.Rel;          % [mM/ms]

J_serca = (c.G.SERCA)*Q10SRCaP^Qpow*Vmax_SRCaP*((yxt(38)/Kmf)^hillSRCaP-(yxt(31)/Kmr)^hillSRCaP)...
    /(1+(yxt(38)/Kmf)^hillSRCaP+(yxt(31)/Kmr)^hillSRCaP); %[44]
J_SRleak = (c.G.Leak)*5.348e-6*(yxt(31)-yxt(36));  % [mM/ms 46]

%% Sodium and Calcium Buffering
ydot(17) = kon_na*yxt(32)*(Bmax_Naj-yxt(17))-koff_na*yxt(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*yxt(33)*(Bmax_Nasl-yxt(18))-koff_na*yxt(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers胞质Ca缓冲液
ydot(19) = kon_tncl*yxt(38)*(Bmax_TnClow-yxt(19))-koff_tncl*yxt(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*yxt(38)*(Bmax_TnChigh-yxt(20)-yxt(21))-koff_tnchca*yxt(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-yxt(20)-yxt(21))-koff_tnchmg*yxt(21);   % TnCHm     [mM/ms]
ydot(22) = kon_cam*yxt(38)*(Bmax_CaM-yxt(22))-koff_cam*yxt(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*yxt(38)*(Bmax_myosin-yxt(23)-yxt(24))-koff_myoca*yxt(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-yxt(23)-yxt(24))-koff_myomg*yxt(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*yxt(38)*(Bmax_SR-yxt(25))-koff_sr*yxt(25);                    % SRB       [mM/ms]
%J_CaB_cytosol = sum(ydot(19:25)); % wrong formulation
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*yxt(36)*(Bmax_SLlowj-yxt(26))-koff_sll*yxt(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*yxt(37)*(Bmax_SLlowsl-yxt(27))-koff_sll*yxt(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*yxt(36)*(Bmax_SLhighj-yxt(28))-koff_slh*yxt(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*yxt(37)*(Bmax_SLhighsl-yxt(29))-koff_slh*yxt(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations SR Ca浓度
ydot(30) = kon_csqn*yxt(31)*(Bmax_Csqn-yxt(30))-koff_csqn*yxt(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations钠浓度
I_Na_tot_junc = I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc+I_NaL_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl+I_NaL_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(yxt(33)-yxt(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(yxt(32)-yxt(33))...
   +J_na_slmyo/Vsl*(yxt(34)-yxt(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) = J_na_slmyo/Vmyo*(yxt(33)-yxt(34));             % [mM/msec] 
% ydot(34)=0;

% Potassium Concentration钾浓度
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur;     % [uA/uF] %SVP: added IKur
% ydot(35) = 0; %-I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) =0; % -I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations钙浓度
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(yxt(37)-yxt(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(yxt(36)-yxt(37))...
    + J_ca_slmyo/Vsl*(yxt(38)-yxt(37))-J_CaB_sl;   % Ca_sl
% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(yxt(37)-yxt(38));
% ydot(38)=0

%% Membrane Potential
I_Na_tot = I_Na_tot_junc+I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk+I_ClCFTR;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;

%用来电流识别
if Id >= 20
    if t >= 100 && mod(t,100) <= 5 % mod 取余数
        Id = 12.5;
    else
        Id = 0;
    end
end

ydot(39) = -(I_tot-Id);

%电流
currents = [I_Na I_to I_kur I_kr I_ks I_ki I_Catot I_ncx I_nak ...
    I_nabk I_cabk I_pca];

% ----- END EC COUPLING MODEL ---------------
if flag
    out=ydot';
else
    
    out=currents;
end
