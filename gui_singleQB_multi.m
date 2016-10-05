function [I tM PLCmax PLCtot DAGmax Ntrpmax Camax Cafreemax]=gui_singleQB_multi(param_set)

C = 2.4;   % conversion factor in molecules/uM
ro = 0.25;  % round onset is at 0.5+ro, so here round(x-ro)=1 if x=0.9. NOT used in stochastic mode.
% eff = 0.9; % collision efficiency

%{
if stoch_flag==0
    global Kcalx;
    global vpkc;
    global Kpkc1_50;
    global Dpip2;
    flagM=0;
    flagG=0;
    flagP=0;
    flagD=0;
    flagT=0;
else

%    Kcalx = 0.4000; % Calcium concentration for 50% activity of pump in mM
%    vpkc = 0.0425; % maximum rate of TRP phosphorylation by PKC (for Apkc = 1)
%    Kpkc1_50 = 145; %DAG amount for 50% DAG-induced activity of PKC on TRP in molecules
%    Dpip2 = 1.1250; % Diffusion constant for PIP2 molecule in membrane in um^2/sec
%    flagM=1;
%    flagG=1;
%    flagP=1;
%    flagD=1;
%    flagT=1;
%end
%}

%{
%flags=[par_flag;flagM;flagG;flagP;flagD;flagT];

%flagM=flags(2);
%flagG=flags(3);
%flagP=flags(4);
%flagD=flags(5);
%flagT=flags(6);


%MAIN FILE: SIMULATES A QUANTUM BUMP

%if strcmp(mode,'single')==1
%    Npoints = 1+time/tstep; 
%end
%if strcmp(mode,'average')==1
    
%end
%}

%{
global flagM;
global flagG;
global flagP;
global flagD;
global flagT;
%}
buffer = load('buffer50.txt');    % nonlinear buffer, calmodulin 0.5mM

%{

##Redundant globals##

%global time; 
%global tstep; 

%global flagM;
%global flagG; 
%global flagP; 
%global flagD; 
%global flagT;

%global tDAGdelay
%global Icalx_sat 
%global Yb_dark;
%global KT;
%global KR;
%global buffer;

%global Dmetarh; global Dgi; global Dga; global Dpip2;
%global alpha1; global alpha2; global alpha3;
%global Dca; global Dmg; global Dna; global Dk; 
%global Smv; 
%global Vmv; global Snk; global Lnk;
%global wca;
%global wmg; global wna; global wk; 
%global Vm;  global F;  
%global Ca_o; global Mg_o; global Na_o; global K_o;   
%global Ca_i; global Mg_i; global Na_i; global K_i; 
%}

 

%% PARAMETERS
%{
%tDAGdelay = variables(1);
global Yb_dark; %= variables(1);
global KT;% = variables(2);
global KR;% = variables(3);
global Icalx_sat;% = variables(4);
global Yb_max;% = variables(5);
global kMA;% = variables(6);
%P1max = variables(8);
global bet;% = variables(7);
%Dgi = variables(10);
%Dga = variables(11);
global Kcalx;% = variables(8);
global vpkc;% = variables(9);
global Kpkc1_50;% = variables(10);
global Dpip2;% = variables(11);

global tDAGdelay; %=10;  %12           %tDAGdelay=tDAGdelay1;  % 10
global P1max; %1.4;     %1.9000;           %P1max=P1max8;
global Dgi; %=1.2;     %1.3750;             %Dgi=Dgi10;
global Dga; %=1.5;     %1.6406;             %Dga=Dga11;

global Mo; % = 1;             % total number of isomerized M molecules
global Gtot; % = 100;         % total number of G protein molecules
global Ptot; % = 100;         % total number of PLC molecules
global PIPo; % = 3e3; %5000   % total number of PIP2 molecules in dark
global Atot; % = 70;          % total number of Arrestin molecules
global Ntot; % = 100;         % total number of NINAC molecules 
global Ntrpo; % = 25;         % total number of TRP channels
global PKCtot; % = 100;       % total number of PKC molecules

global Kninac_max; % = 3;     % maximum value for association constant NINAC-Arr2 in uM^(-1)
global Kcam; % = 0.01;  % 0.05  % conc of Ca for 50% activity of calmodulin in mM
global beta1; % = 25;         % exponential parameter for the reaction rate of NINAC + Arr -> NINAC*Arr
                    % higher kMA gives longer tau_M liftime, smaller Kcam
                    % makes it more sensitive to Ca influx - the aim is to
                    % have tau_M~25ms, but latency time ~45ms
global tGDP; % = 5;     % 4      % time for GDP-GTP exchange in ms                                                      
global t1; % = 1;     % 4        % delay time for Galpha-GTP protein activation in ms  
global tpi; % = 1;  % (9) 8   % time for reaction GPLC+PIP2->DAG in dark in ms
global Kpi_50; % = 0.05; % conc of Ca for 50% activity of PLC in mM
global beta2; % = -2; %  1.3 % exp param for rate of PLC+PIP2->PLC+DAG
global tPdark; % = 100;  %100 % time constant for decay of GPLC in dark in ms
global tDdark; % = 150;       % time constant for decay of DAG in dark in ms
global beta3; % = 3.5; % 3.5   % exp param for reaction rate of GPLC*+GAP->G+PLC
global beta4; % = 4.5; % 4.5  % exp param for rate of DAG + DGK -> PA ... -> PIP2
global K50_gap; % = 0.1;      % concentration of calcium for 50% activity of GAP in mM
global K50_dgk; % = 0.5;      % concentration of calcium for 50% activity of DGK in mM
global t2; % = 20;            % delay time for Ca-induced activation of GAP
global t3; % = 30;            % delay time for Ca-induced activation of DGK
global t4; % = 40;            % delay time for action of PKC on TRP
global vph; % = 0.4;          % rate of TRP relaxation in dark (dephosphorylation) in ms^(-1)
global beta5; % = 6;      % exponential parameter for the reaction TRPp -> TRP

global Dmetarh; % = 0.0; % Diffusion const of active metarhodopsin um^2/s

global alpha1; % = 1.960e-3;
global alpha2; % = 2.005e-3;
global alpha3; % = 1.706e-3;
 alpha1 = 4*3.1415/(log(4*(Dmetarh+Dgi)*tcharac1/rho1^2)-1.15);
 alpha2 = 4*3.1415/(log(4*Dga*tcharac2/rho2^2)-1.15);
 alpha3 = 4*3.1415/(log(4*Dpip2*tcharac3/rho3^2)-1.15);

global Dca; % = 220;      % Diffusion constant for Ca ion in microvillus in um^2/sec
global Dmg; % = 200;      % Diffusion constant for Mg ion in microvillus in um^2/sec
global Dna; % = 650;      % Diffusion constant for Na ion in microvillus in um^2/sec
global Dk; % = 1000;      % Diffusion constant for K ion in microvillus in um^2/sec

global Smv; % = 0.27;         % Surface of microvillus in um^2
global Vmv; % = 4.24e-18;     % Volume of microvillus in L
global Snk; % = 9.6e-4;       % Surface of microvillus neck in um^2
global Lnk; % = 0.06;         % Length of microvillus neck in um


global n; % = 4;         % Number of subunits in TRP channel (MWC model parameter)
global Kcamtrp; % = 0.2; %
global Kpkc2_50; % = 1.0; % Calcium conc for 50% Ca-induced activity of PKC on TRP in mM

global wca; % = 0.877; %0.9*0.88+0.1*0.85; Fraction of total Permeability due to Ca
global wmg; % = 0.101; %0.9*0.10+0.1*0.11; Fraction of total Permeability due to Mg
global wna; % = 0.011; %0.9*0.01+0.1*0.02; Fraction of total Permeability due to Na
global wk; %  = 0.011; %0.9*0.01+0.1*0.02; Fraction of total Permeability due to K

global Vm; % = -0.070;        % Holding membrane potential in V
global F; % = 96500;          % Faraday constant in C/mol

global Ca_o;% = 1.5;     % Physiological Ca extracellular concentration in mM
global Mg_o; % = 4;           % Physiological Mg extracellular concentration in mM
global Na_o; % = 120;         % Physiological Na extracellular concentration in mM
global K_o; % = 5.0;          % Physiological K extracellular concentration in mM

global Ca_i; % = 1.6e-4;% 16;  % Physiological Ca intracellular concentration in mM
global Mg_i; % = 3;        % Physiological Mg intracellular concentration in mM
global Na_i; % = 8.0;      % Physiological Na intracellular concentration in mM
global K_i; % = 140;       % Physiological K intracellular concentration in mM

%{
##Redundant parameters##

 tcharac1=0.01;  % s
 tcharac2=0.01;  % s
 tcharac3=0.02;  % s
 d_M = 0.003; % um
 d_G = 0.007; % um
 d_Ga = 0.005;
 d_PLC = 0.007;
 d_pip2 = 0.001;
 rho1 = (d_M+d_G)/2;
 rho2 = (d_Ga+d_PLC)/2;
 rho3 = (d_PLC+d_pip2)/2;
 alpha1 = 1e-3*4*3.1415/(log(4*(Dmetarh+Dgi)*tcharac1/rho1^2)-1.15);
 alpha2 = 1e-3*4*3.1415/(log(4*Dga*tcharac2/rho2^2)-1.15);
 alpha3 = 1e-3*4*3.1415/(log(4*Dpip2*tcharac3/rho3^2)-1.15);

 kMA = 4e-3; 9, 2.5 for 200ms, reaction rate for binding M-Arr2 in (molecule*ms)^(-1)

 Dgi = 1.2; %1.2 % Diffusion constant for inactive G protein in membrane in um^2/sec
 Dga = 1.5;      % Diffusion constant for active Galpha protein in membrane in um^2/sec
 Dplc = 0; % part of INAD-complex, practically not moving
 Dgplc = 0; % same 

 tDAGdelay = 1; % ms DAG delay to bind to TRP
 Yb_dark = 0.3e-6; % 8e-7 Basal activity for allosteric transition in dark (MWC model parameter)
 Yb_max = 3e-6;  % 8e-6 Maximum CaM/Ca-induced basal activity for allosteric transition 
 KT = 0.005; % 5e-3; % Affinity to ligand for T-state (MWC model parameter) in molecules^(-1)
 KR = 0.285; % 0.56   % Affinity to ligand for R-state (MWC model parameter) in molecules^(-1)
 Icalx_sat = 5.0; %1.6   % Saturation current for exchange CalX pump in pA

 P1max = 1.40;  % (1.8) 1.39     % Maximum permeability of a single TRP channel in um/s

%}
%}

for k=1:77
    eval(param_set{k});
end

p = nextpow2(time/tstep);
Npoints = 2^p;

tup = (Npoints-1)*tstep;       % ms

t = 0:tstep:tup;               % ms
s = length(t);

itDAGdl = round(tDAGdelay/tstep);

F=96500;

%% FUNCTIONS AND VARIABLES
%{ 
 tcharac1=0.01;  % s
 tcharac2=0.01;  % s
 tcharac3=0.01;  % s
 d_M = 0.003; % um
 d_G = 0.007; % um
 d_Ga = 0.005;
 d_PLC = 0.007;
 d_pip2 = 0.001;
 rho1 = (d_M+d_G)/2;
 rho2 = (d_Ga+d_PLC)/2;
 rho3 = (d_PLC+d_pip2)/2;
 alpha1 = 4*3.1415/(log(4*Dgi*tcharac1/rho1^2)-1.15)
 alpha2 = 4*3.1415/(log(4*Dga*tcharac2/rho2^2)-1.15)
 alpha3 = 4*3.1415/(log(4*Dpip2*tcharac3/rho3^2)-1.15)
%}

%buffer = load('buffer50.txt');    % nonlinear buffer, calmodulin 0.5mM

%ddt_1 = exp(-t/t1)/t1;            % delay stage functions
%d2 = 1-exp(-t/t2);
%d3 = 1-exp(-t/t3);
%d4 = 1-exp(-t/t4);

% M*
M = zeros(1,s);         M(1) = Mo;
MM = M;       
Arr_free = zeros(1,s);
tau_M  = zeros(1,s);
Arr_free(1) = 1;
tau_M(1) = 1/(Arr_free(1)*kMA);

% G*
Gc1 = zeros(1,s);                    % "c" label: continous
Gc2 = zeros(1,s);
Gd = zeros(1,s);                     % "d" label: discrete

% GPLC*
GPLCc1 = zeros(1,s);
GPLCc2 = zeros(1,s);
GPLCd = zeros(1,s);
PLCtot = 0;

% DAG
PIP2loss = zeros(1,s);
DAGc1 = zeros(1,s);
DAGc2 = zeros(1,s);
DAGd = zeros(1,s);
dDAGdt = zeros(1,s);

% TRP
CH = zeros(1,Ntrpo);    % array of channels (0 or 1) for stochastic mode
Ntrp = zeros(1,s);
Nact = zeros(1,s);      Nact(1) = Ntrpo;
Nactc1 = zeros(1,s); 
Nactc2 = zeros(1,s);
%Yb = zeros(1,s);       
Yb_old = Yb_dark;   

% PLC, GAP, DGK, PKC
Aplc = zeros(1,s);
Agap = zeros(1,s);  
Adgk = zeros(1,s);  
Apkc = zeros(1,s);
% activation functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%activation1 = zeros(1,s);
%activation2 = zeros(1,s);
%activation3 = zeros(1,s);
%activation4 = zeros(1,s);
%activation5 = zeros(1,s);
%Knina= zeros(1,s);
% ...........................................................
% Ca2+, Mg2+, Na+, K+
Ca = zeros(1,s);        Ca(1) = Ca_i;   
dCadt = zeros(1,s);
Cafree = zeros(1,s);    Cafree(1) = Ca_i;
dCafreedt = zeros(1,s);

Mg = zeros(1,s);        Mg(1) = Mg_i;
Na = zeros(1,s);        Na(1) = Na_i;
K = zeros(1,s);         K(1) = K_i;

% Itrp, Icalx, Iq
Itrp = zeros(1,s);
Icalx = zeros(1,s);
Ica = zeros(1,s);
Img = zeros(1,s);
Ina = zeros(1,s);
Ik = zeros(1,s);

% Threshold
T = zeros(1,s);   % dynamic threshold (= DAG amount needed for 1 open channel)

%istim2 = 0; 
%istim3 = 0;   % times used to reset convolutions in case of value coming back to zero.
%istim4 = 0;

%---------------------------------
Xg = zeros(1,s);
Xgap = zeros(1,s);
Xdgk = zeros(1,s);
Xpkc = zeros(1,s);

vcoll_PIP2=zeros(1,s);
tau_PIP2=zeros(1,s);
vPIP2=zeros(1,s);
%----------------------------------

for i=2:s
    
%RHODOPSIN DYNAMICS
Acam = Ca(i-1)/(Ca(i-1)+Kcam);

Kninac = kninac_max*exp(-beta1*Acam);  
%Kninac = kninac_max*(1-Acam); 
%Knina(i)=Kninac;
%activation1(i)=beta1*Acam; %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a = Kninac;
b = 1+ Kninac*(Ntot-Atot)/C;
c = - Atot/C;
Afree = round(C*(-b+sqrt(b^2-4*a*c))/(2*a)-ro);
Arr_free(i) = Afree;
tau_M(i) = 1/(kMA*Afree);

if flagM==0
%  dMdt = - kMA*Afree*M(i-1); 
%  M(i) = M(i-1) + dMdt*tstep;
% 1 new stuff
    dMdt = - kMA*Afree*MM(i-1);
    MM(i) = MM(i-1) + dMdt*tstep;
    M(i) = round(MM(i)-ro);
%====================================
else
    deltaM = - poisson(kMA*Afree*M(i-1)*tstep);
    M(i) = M(i-1) + deltaM;
end

% multiple photon option: to simulate another photon at time t_photon
%{
t_photon = 200;
if i == t_photon/tstep
    
    M(i) = Mo;
    
end
%}  

    
% reset test
if M(i)<1e-4
    M(i)=0;
end

%CASCADE DYNAMICS
% DYNAMICS

%c1 = convolution(M, ddt_1, 0, i-1, tstep);

%c2 = convolution(dCadt.*(K50_gap./(K50_gap+Ca).^2),d2,istim2,i-1,tstep);
%c3 = convolution(dCadt.*(K50_dgk./(K50_dgk+Ca).^2),d3,istim3,i-1,tstep);

vcoll_Gi = alpha1*((Dmetarh+Dgi)/Smv)*(Gtot - Gd(i-1) - GPLCd(i-1));        % rate of G-M* coll ms^(-1)
vcoll_Ga = alpha2*(Dga/Smv)*(Ptot-GPLCd(i-1))/((1+sqrt(GPLCd(i-1)/pi))^2);  % rate of G*-PLC coll ms^(-1)
vcoll_PI = alpha3*(Dpip2/Smv)*(PIPo - PIP2loss(i-1)); % -DAGd(i-1));  % rate of PIP2-GPLC* coll ms^(-1)
vcoll_PIP2(i) = vcoll_PI;

Aplc(i) = Cafree(i-1)/(Cafree(i-1)+Kpi_50);
tPI = tpi*exp(beta3*Aplc(i)); % 
tau_PIP2(i) = tPI;
%activation3(i)=beta3*Aplc(i);  %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

vGi = vcoll_Gi/(1+vcoll_Gi*tGDP); % rate of Galpha generation in ms^(-1)
vGa = vcoll_Ga;                   % rate of GPLC generation in ms^(-1)
vPI = vcoll_PI/(1+vcoll_PI*tPI);  % rate of DAG generation in ms^(-1)
vPIP2(i) = vPI;

%Agap(i) = c2(i-1);
%Adgk(i) = c3(i-1);

%----------------------------------------

Agap(i) = Agap(i-1) + Xgap(i-1)*tstep;
Sgap = dCadt(i-1)*(K50_gap/(K50_gap+Ca(i-1))^2);
Xgap(i) = Xgap(i-1) + tstep*(Sgap - Xgap(i-1))/t2;

Adgk(i) = Adgk(i-1) + Xdgk(i-1)*tstep;
Sdgk = dCadt(i-1)*(K50_dgk/(K50_dgk+Ca(i-1))^2);
Xdgk(i) = Xdgk(i-1) + tstep*(Sdgk - Xdgk(i-1))/t3;


%----------------------------------------

%consistency test
if Agap(i)<0
    Agap(i)=0;
    %istim2 = i;
end
if Adgk(i)<0
    Adgk(i)=0;
    %istim3 = i;
end

tP = tPdark*exp(-beta2*Agap(i));
tD = tDdark*exp(-beta4*Adgk(i));
%activation2(i) = beta2*Agap(i); %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%activation4(i) = beta4*Agap(i);

if flagG == 0
    %----------------------------------
    Xg(i) = Xg(i-1) + tstep*(vGi*M(i-1)-Xg(i-1))/t1;
    dGdt1 = Xg(i);
    %----------------------------------
    %dGdt1 = vGi*c1(i-1);
    dGdt2 = vGa*Gd(i-1);
    Gc1(i) = Gc1(i-1) + dGdt1*tstep;
    Gc2(i) = Gc2(i-1) + dGdt2*tstep;
%2 quantised or continuous deterministic values (G) ?????????????????
    Gd(i) = round(Gc1(i)-ro) - round(Gc2(i)-ro);
%   Gd(i) = Gc1(i) - Gc2(i);
else
    %----------------------------------
    Xg(i) = Xg(i-1) + tstep*(vGi*M(i-1)-Xg(i-1))/t1;
    dGdt1 = Xg(i);
    %----------------------------------
    %dGdt1 = vGi*c1(i-1);
    dGdt2 = vGa*Gd(i-1);  
    deltaG1 = poisson(dGdt1*tstep);
    deltaG2 = poisson(dGdt2*tstep);
    if deltaG2 > Gd(i-1)+deltaG1
        deltaG2 = Gd(i-1)+deltaG1;
    end
    Gd(i) = Gd(i-1) + deltaG1 - deltaG2;
end

%consistency test
if Gd(i)<0
    Gd(i)=0;
end

%reset test
if Gd(i)==0 && M(i)==0
    Gc1(i) = 0;
    Gc2(i) = 0;
end

if flagP == 0
   dPdt1 = dGdt2;
   dPdt2 = (1/tP)*GPLCd(i-1); 
   GPLCc1(i) = GPLCc1(i-1) + dPdt1*tstep;
   GPLCc2(i) = GPLCc2(i-1) + dPdt2*tstep;
% 3 quantised or continuous deterministic values (G-PLC) ?????????????????
   GPLCd(i) = round(GPLCc1(i)-ro) - round(GPLCc2(i)-ro); 
%   GPLCd(i) = GPLCc1(i) - GPLCc2(i); 
% new line for counting total active PLC
  PLCtot=PLCtot+(round(GPLCc1(i)-ro)-round(GPLCc1(i-1)-ro));
% ======================================
else
   dPdt2 = (1/tP)*GPLCd(i-1); 
   deltaP1 = deltaG2;  
   deltaP2 = poisson(dPdt2*tstep);
   GPLCd(i) = GPLCd(i-1) + deltaP1 - deltaP2; 
% new line for counting total active PLC
  PLCtot=PLCtot+deltaP1;
% ======================================
end

%consistency test
if GPLCd(i)<0
    GPLCd(i)=0;
end

%reset test
if GPLCd(i)==0 && M(i)==0
    GPLCc1(i) = 0;
    GPLCc2(i) = 0;
end

if flagD == 0
    dDdt1 = vPI*GPLCd(i-1);
    dDdt2 = (1/tD)*DAGd(i-1);
    DAGc1(i) = DAGc1(i-1) + dDdt1*tstep;
    DAGc2(i) = DAGc2(i-1) + dDdt2*tstep;
%old    DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i-1)-ro);
% 4 quantised or continuous deterministic values (DAG) ?????????????????
    DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i)-ro);
%    DAGd(i) = DAGc1(i) - DAGc2(i);
    dDAGdt(i) = (DAGd(i)-DAGd(i-1))/tstep;
else
    dDdt1 = vPI*GPLCd(i-1);
    dDdt2 = (1/tD)*DAGd(i-1);
    DAGc1(i) = DAGc1(i-1) + poisson(dDdt1*tstep);
    DAGc2(i) = DAGc2(i-1) + poisson(dDdt2*tstep);
%    DAGd(i) = round(DAGc1(i)-ro) - round(DAGc2(i-1)-ro);
    DAGd(i) = DAGc1(i) - DAGc2(i);
end
% new line for counting consumed PIP2
%    PIP2loss(i)=PIP2loss(i-1)+DAGd(i);
     PIP2loss(i)=DAGc1(i);
% ======================================

%consistency test

if DAGd(i)<0
    DAGd(i)=0;
end

%reset test
if DAGd(i)==0 && M(i)==0
    DAGc1(i) = 0;
    DAGc2(i) = 0;
end


%CHANNEL DYNAMICS
% new ch dynamics
L=0;
if i>itDAGdl
    L=DAGd(i-itDAGdl);
end
% =========================
%L = DAGd(i-1);
Yo = Yb_old;

Atrp = ((1+KR*L).^n )./( (1+KR*L).^n + (1/Yo)*(1+KT*L).^n );

%bet = 4;    %  2     % mean rate of channel closing (only for stochastic mode)
alf = bet*Atrp/(1-Atrp);   % rate of channel opening

if flagT == 0
% 5 quantised or continuous deterministic values (Ntrp) ?????????????????
    Ntrp(i) = round(Nact(i-1)*Atrp-ro); 
%    Ntrp(i) = Nact(i-1)*Atrp; 
else
    alfatau = alf*tstep;
    betatau = bet*tstep;
    for j=1:Nact(i-1)
        if CH(j)==0
%            CH(j) = poisson(alf*tstep);
%            if CH(j)>0  % 1
%                CH(j) = 1;
%            end   
            if alfatau > rand(1)  % 1
                CH(j) = 1;
            end
        else 
%            CH(j) = 1 - poisson(bet*tstep);
%            if CH(j)<1  % 0
%                CH(j)=0;
%            end
            if betatau > rand(1)  % 1
                CH(j) = 0;
            end
        end
    end  
      
    for j=Nact(i-1)+1:Ntrpo
        CH(j)=0;
    end
    
    Ntrp(i) = sum(CH);
    
end

%if Ntrp(i)>0.5
        %check_dis2=Ntrp(i)
%end
    
Acam = Ca(i-1)/(Ca(i-1)+Kcam);
%new function for +ve feedback via TRP activation 
Acamtrp = Ca(i-1)/(Ca(i-1)+Kcamtrp); % new line
Yb_new = Yb_dark + (Yb_max-Yb_dark)*Acamtrp; %Acam

vrelax = vph*exp(-beta5*Acam);
%activation5(i) = beta5*Agap(i); %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if Apkc(i-1)>0
dNdt1 = vpkc*PKCtot*Apkc(i-1)*Nact(i-1);  %  Ntrp(i-1)
else
dNdt1=0;
end
dNdt2 = vrelax*(Ntrpo - Nact(i-1));

Nactc1(i) = Nactc1(i-1) + dNdt1*tstep;
Nactc2(i) = Nactc2(i-1) + dNdt2*tstep;
%Nactc1(i) = Nactc1(i-1) + poisson(dNdt1*tstep);
%Nactc2(i) = Nactc2(i-1) + poisson(dNdt2*tstep);

if flagT == 0
% 6 
    Nact(i) = Ntrpo - round(Nactc1(i)-ro) + round(Nactc2(i)-ro);
%    Nact(i) = Ntrpo - Nactc1(i) + Nactc2(i);
else
    Nact(i) = Ntrpo - round(Nactc1(i)-ro) + round(Nactc2(i)-ro);
end

%A1 = (DAGd./(Kpkc1_50 + DAGd));
%A2 = (Cafree./(Kpkc2_50 + Cafree));
%dApkcdt = A1.*dCafreedt*Kpkc2_50./(Kpkc2_50+Cafree).^2 + A2.*dDAGdt*Kpkc1_50./(Kpkc1_50+DAGd).^2 ;
%c4 = convolution(dApkcdt,d4,istim4,i-1,tstep);

%Apkc(i) = c4(i-1);

%--------------------------------------------

Apkc(i) = Apkc(i-1) + Xpkc(i-1)*tstep;
A1 = (DAGd(i-1)/(Kpkc1_50 + DAGd(i-1)));
A2 = (Cafree(i-1)/(Kpkc2_50 + Cafree(i-1)));
Spkc = A1*dCafreedt(i-1)*Kpkc2_50/(Kpkc2_50+Cafree(i-1))^2 + A2*dDAGdt(i-1)*Kpkc1_50/(Kpkc1_50+DAGd(i-1))^2 ;
Xpkc(i) = Xpkc(i-1) + tstep*(Spkc - Xpkc(i-1))/t4;
%--------------------------------------------

%consistency test
if Apkc(i)<0 
    Apkc(i)=0;
    %istim4 = i; 
end


%threshold update
alpha = (Nact(i)*Yb_new)^(1/4);
if alpha==0
    T(i)=inf;
else
    T(i) = (1-alpha)/(alpha*KR - KT);
end



%CURRENT DYNAMICS

Icalx(i) = -Icalx_sat*Cafree(i-1)/(Cafree(i-1) + Kcalx);

Ica(i) = ghk_current(Ntrp(i-1),Cafree(i-1),Ca_o,P1max,wca,2,Vm,Smv);
Img(i) = ghk_current(Ntrp(i-1),Mg(i-1),Mg_o,P1max,wmg,2,Vm,Smv);
Ina(i) = ghk_current(Ntrp(i-1),Na(i-1),Na_o,P1max,wna,1,Vm,Smv);
Ik(i) = ghk_current(Ntrp(i-1),K(i-1),K_o,P1max,wk,1,Vm,Smv);

Itrp(i) = Ica(i) + Img(i) + Ina(i) + Ik(i);

Ica_leak = leak_current(2, Dca, Cafree(i-1), Ca_i, Lnk, Snk);
Img_leak = leak_current(2, Dmg, Mg(i-1), Mg_i, Lnk, Snk);
Ina_leak = leak_current(1, Dna, Na(i-1), Na_i, Lnk, Snk);
Ik_leak = leak_current(1, Dk, K(i-1), K_i, Lnk, Snk);


%ION DYNAMICS

dCadt(i) = 1e-12*(-Ica(i-1)/2 + Icalx(i-1) - Ica_leak/2)/(F*Vmv);      % in mmol/(L*msec)
Ca(i) = Ca(i-1) + dCadt(i)*tstep;                                      % in mM
  if Ca(i)<0
      Ca(i)=0;
  end
i_buff = uint16(round(1e3*Ca(i)) + 1);                     
B = 1+buffer(i_buff); %/2;    % buffer power for 0.5 mM (0.25mM) Calmodulin)
Cafree(i) = Ca(i)/B;          % in mM
dCafreedt(i) = (Cafree(i)-Cafree(i-1))/tstep;

dMgdt = 1e-12*(-Img(i-1) - Img_leak)/(2*F*Vmv);                     
Mg(i) = Mg(i-1) + dMgdt*tstep;              

dNadt = 1e-12*(-Ina(i-1) - 3*Icalx(i-1) - Ina_leak)/(F*Vmv);    
Na(i) = Na(i-1) + dNadt*tstep;  

dKdt = 1e-12*(-Ik(i-1) - Ik_leak)/(F*Vmv);                  
K(i) = K(i-1) + dKdt*tstep;     

%consistency test
if Ca(i)<0
    Ca(i)=0;
    Cafree(i)=0;
end

Yb_old=Yb_new;

end


%% OUTPUT
%noise=0.20; % 0.2
%Iclean = Itrp+Icalx;
%in = length(Iclean);
%Inoise = noise*randn(1,in);
%I = Iclean + Inoise;
I = Itrp+Icalx;
%figure(50)
%plot(t,Itrp,'k',t,Icalx,'b',t,Ica_leak,'r')

tM = sum(M)*tstep;
PLCmax = max(GPLCd);
Ntrpmax = max(Ntrp);
Camax = max(Ca);
Cafreemax = max(Cafree);
DAGmax = max(DAGd);

Imax = max(-I); %#ok<NASGU>

%PIP2lossout = PIP2loss(end);

%% PLOT
% % plot the current for just one single run =============   
% % (1) filtering function ===============================
%     fcut = 100;
%     N = Npoints;
%     kspan = 1+uint16(1000/(3.1415*fcut*tstep));
%     If(1:kspan-1) = 0;
%     iss = 1:kspan-1;
%     ksfloat = single(kspan);
%     aa=1-single(iss)./ksfloat;
%     cc=2/ksfloat;
%          for ns = kspan:N
%               sumI = I(ns);
%                  for is = 1:kspan-1
%                    sumI = sumI + aa(is)*I(ns-is);
%                  end
%               If(ns)=cc*sumI;
%           end
%     figure(6)
%     hold on
%     plot(t,I,'b',t,If,'r')
%========================================================
%{
bump = load('bump_data0.txt');
QQ = bump(:,1);
num2 = length(QQ); % 65;
TT = 0:1:(num2-1);
T = 0:tstep:(num2-1);
istep = 1/tstep;

for i2 = 1:num2
    ir = (i2-1)*istep + 1;
    result(i2) = I(ir);  % reduce the numbe of points to exp
    res_Itrp(i2) = Itrp(ir);
    res_Icalx(i2) = Icalx(ir);
    res_Icaleak(i2) = Icaleak(ir)
end

hold on
plot(TT,QQ,'or',T,I,'k',TT,result,'x', T,Itrp,'g',T,Icalx,'b',T,Icaleak,'r')
hold off
%}

%figure(6)
%%subplot(3,1,1),plot(t,activation1, subplot(3,1,2), plot(t,Knina);subplot(3,1,3)
%plot(t,Arr_free,'k',t,tau_M,'r');
%xlim([0 120])

%{
Plot the activation functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure(6)
subplot(5,1,1), plot(t,activation1);
title('activity1');
subplot(5,1,2), plot(t,activation2);
subplot(5,1,3), plot(t,activation3);
subplot(5,1,4), plot(t,activation4);
subplot(5,1,5), plot(t,activation5);
%xlim([0 xl])
%}

%{
% print for a determ single run
hold on
figure(6)
xl = 200;
%plot(handles.axes4,TT,QQ,'or',T,Is,'k',TT,result,'x');
subplot(6,1,1), plot(t,M,'k');
title('Metarhodopsin activity');
xlim([0 xl])
ylim([0 1.1])

subplot(6,1,2), plot(t,Gd,'r',t,GPLCd,'k');
title('G* and GPLC*');
xlim([0 xl])

subplot(6,1,3), plot(t,PIP2loss,'r',t,DAGd,'k');
title('PIP2loss, DAG(t)');
xlim([0 xl])
ylim([0 150])

subplot(6,1,4), plot(t,Ntrp,'r',t,Nact,'b');
title('TRP channels active (blue) & open (red)');
xlim([0 xl])
ylim([0 27])

subplot(6,1,5), plot(t,I,'b'); %,t,Ifilter,'b')
title('Filtered Photocurrent');
xlim([0 xl])
ylim([-10 0.2])

subplot(6,1,6), plot(t,Ca,'r',t,Cafree,'k');
title('Calcium concentration (mM)');
%xlim([0 100])
xlim([0 xl])
ylim([0 2.4])
hold off
%}


