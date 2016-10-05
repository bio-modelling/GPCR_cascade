 function gui_averageQBall(param_set, m_set, procflag)
 
starttime=clock;
 
 for i=1:4
     eval(m_set{i});
 end

eval(param_set{76});
eval(param_set{77});
eval(param_set{78});

%{
%global Ilat 
%global Ibump 
%global tmin 
%global fcut 
%global N
%global tstep; 
%global T12

%global av_tM;
%global av_peakPLC;
%global av_PLCtot;
%global av_peakDAG;
%global av_peakNtrp;
%global av_peakCa;

%global buffer;
%}
%{
global Dmetarh; global Dgi; global Dga; global Dpip2;
global alpha1; global alpha2; global alpha3;
global Dca; global Dmg; global Dna; global Dk; 
global Smv; global Vmv; global Snk; global Lnk;

global wca; global wmg; global wna; global wk; 
global Vm;  global F;  
global Ca_o; global Mg_o; global Na_o; global K_o;   
global Ca_i; global Mg_i; global Na_i; global K_i; 

%tstep = 0.1;  % ms
runs = r;

Dmetarh = 0.0; % Diffusion const of active metarhodopsin um^2/s
Dgi = 1.2; %1.2 % Diffusion constant for inactive G protein in membrane in um^2/sec
Dga = 1.5;      % Diffusion constant for active Galpha protein in membrane in um^2/sec
% Dplc = 0; % part of INAD-complex, practically not moving
% Dgplc = 0; % same
Dpip2 = 1.0; % 1  % Diffusion constant for PIP2 molecule in membrane in um^2/sec
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
 
Dca = 220;      % Diffusion constant for Ca ion in microvillus in um^2/sec
Dmg = 200;      % Diffusion constant for Mg ion in microvillus in um^2/sec
Dna = 650;      % Diffusion constant for Na ion in microvillus in um^2/sec
Dk = 1000;      % Diffusion constant for K ion in microvillus in um^2/sec

Smv = 0.27;         % Surface of microvillus in um^2
Vmv = 4.24e-18;     % Volume of microvillus in L
Snk = 9.6e-4;       % Surface of microvillus neck in um^2
Lnk = 0.06;         % Length of microvillus neck in um

wca = 0.9*0.88+0.1*0.85;   % Fraction of total Permeability due to Calcium
wmg = 0.9*0.10+0.1*0.11;   % Fraction of total Permeability due to Magnesium
wna = 0.9*0.01+0.1*0.02;   % Fraction of total Permeability due to Sodium
wk = 0.9*0.01+0.1*0.02;    % Fraction of total Permeability due to Potassium

Vm = -0.070;        % Holding membrane potential in V
F = 96500;          % Faraday constant in C/mol

Ca_o = 1.5;     % Physiological Ca extracellular concentration in mM
Mg_o = 4;           % Physiological Mg extracellular concentration in mM
Na_o = 120;         % Physiological Na extracellular concentration in mM
K_o = 5.0;          % Physiological K extracellular concentration in mM

Ca_i = 0.00016;% 16;  % Physiological Ca intracellular concentration in mM
Mg_i = 3;        % Physiological Mg intracellular concentration in mM
Na_i = 8.0;      % Physiological Na intracellular concentration in mM
K_i = 140;       % Physiological K intracellular concentration in mM
%}

% nonlinear buffer, calmodulin 0.5mM
%buffer = load('buffer50.txt');    
%runs = r;

ip = nextpow2(time/tstep);
N = 2^ip; 
tup = (N-1)*tstep;       % ms
t = 0:tstep:tup;

Iset = zeros(1,N);   % original noisy bumps
Iset1 = zeros(1,N);  % filtered bumps 
Iset2 = zeros(1,N);  % filtered and shifted bumps 
Imaxset = 0; 
L = 0;
t12 = 0;
time1 = 0;
time2 = 0;
time3 = 0;
av_tM = 0;
av_peakPLC = 0;
av_PLCtot = 0;
av_peakDAG = 0;
av_peakNtrp = 0;
av_peakCa = 0;
av_peakCafree = 0;

w = 1; p = 1; s = 1; q = 1;
neff = runs;            % number of bumps "not failed" (Q.E. = neff/runs)

%% Run and filter
    
%{    
if runs==1
    
    for m=1:runs
        [I tM PLCmax PLCtot DAGmax Ntrpmax Camax Cafreemax]=gui_singleQB_multi(param_set); 
        Iset(m,:)=I;
        tMset(m,:)=tM;
        PLCmaxset(m,:)=PLCmax;
        DAGmaxset(m,:)=DAGmax;
        Ntrpmaxset(m,:)=Ntrpmax;
        Camaxset(m,:)=Camax;
        Cafreemaxset(m,:)=Cafreemax;        
    end
else
%}
if procflag==2
    matlabpool open local 2
end

if procflag==4
    matlabpool open local 4
end
    parfor m=1:runs
[I tM PLCmax PLCtot DAGmax Ntrpmax Camax Cafreemax]=gui_singleQB_multi(param_set); 

           Iset(m,:)=I;
           tMset(m,:)=tM;
           PLCmaxset(m)=PLCmax;
           PLCtotset(m)=PLCtot;
           DAGmaxset(m)=DAGmax;
           Ntrpmaxset(m)=Ntrpmax;
           Camaxset(m)=Camax;
           Cafreemaxset(m)=Cafreemax;
    end
%end
if procflag==2||procflag==4
matlabpool close
end

    kspan = 1+uint16(1000/(3.1415*fcut*tstep));
    If(1:kspan-1) = 0;
    iss = 1:kspan-1;
    ksfloat = single(kspan);
    aa=1-single(iss)./ksfloat;
    cc=2/ksfloat;

for m=1:runs
    Inew = Iset(m,:);
% (1) filtering function ===============================
%       If = Filter(Inew,t,fcut); 
          for ns = kspan:N
              sumI = Inew(ns);
                 for is = 1:kspan-1
                   sumI = sumI + aa(is)*Inew(ns-is);
                 end
              If(ns)=cc*sumI;
          end
%}       
% one off =========================================
%figure(7)
%plot(t,If,'b'); %,t,I,'b')
%title('Filtered Photocurrent');
%======================================

   [Imax,k50] = max(-If);
%   Imax %#ok<NOPRT>
   Imaxall(m) = Imax;

 % select the QBs only with Imax greater than some min value above the
 % noise, set by Ibump
   if Imax > Ibump
       lat = 0;
       t1 = 0;
       t2 = 0;
       for k=1:N
           if -If(k) > Ilat && lat == 0 && k > tmin/tstep
               lat =(k-1)*tstep;
           end
           
           if -If(k) >= Imax/2 && t1 == 0
               t1 = (k-1)*tstep;
           end
           
           if -If(k) <= Imax/2 && k>k50 && t2 == 0
               t2 = (k-1)*tstep;
           end
       end
       
       Iset1(w,:) = If;
       Imaxset(p) = Imax;  %#ok<*AGROW>
       L(s) = lat; 
       time1(s) = t1-lat;
       tpeak = (k50-1)*tstep;
       time2(s) = tpeak - t1;
       time3(s) = t2-tpeak;
       t12(q) = (t1+t2)/2; 
       av_tM = av_tM + tMset(w);
       av_peakPLC = av_peakPLC + PLCmaxset(w);
       av_PLCtot = av_PLCtot + PLCtotset(w);
       av_peakDAG = av_peakDAG + DAGmaxset(w);
       av_peakNtrp = av_peakNtrp + Ntrpmaxset(w);
       av_peakCa = av_peakCa + Camaxset(w);
       av_peakCafree = av_peakCafree + Cafreemaxset(w);
       w = w+1; p = p+1; s = s+1; q = q+1;
   else
  
       neff = neff-1; 
   %    av_tM = av_tM + tM;
       
   end
   
     
end


%{    
%parfor m=1:runs
    
%   m    %#ok<NOPRT>
[t I]=gui_singleQB(param_set);    
%[I tM PLCmax PLCtot DAGmax Ntrpmax Camax Cafreemax t] = par_singleQB_all(time,tstep,variables,flags);
   
   Iset(m,:)=I;
   If = Filter(I,t,fcut);   
   
   [Imax,k50] = max(-If);
%   Imax %#ok<NOPRT>
   Imaxall(m) = Imax;
   
   if Imax > Ibump
       lat = 0;
       t1 = 0;
       t2 = 0;
       for k=1:N
           if -If(k) > Ilat && lat == 0 && k > tmin/tstep
               lat =(k-1)*tstep;
           end
           
           if -If(k) >= Imax/2 && t1 == 0
               t1 = (k-1)*tstep;
           end
           
           if -If(k) <= Imax/2 && k>k50 && t2 == 0
               t2 = (k-1)*tstep;
           end
       end
       
       Iset1(w,:) = If;
       Imaxset(p) = Imax;  %#ok<*AGROW>
       L(s) = lat; 
       time1(s) = t1-lat;
       tpeak = (k50-1)*tstep;
       time2(s) = tpeak - t1;
       time3(s) = t2-tpeak;
       t12(q) = (t1+t2)/2; 
       av_tM = av_tM + tM;
       av_peakPLC = av_peakPLC + PLCmax;
       av_PLCtot = av_PLCtot + PLCtot;
       av_peakDAG = av_peakDAG + DAGmax;
       av_peakNtrp = av_peakNtrp + Ntrpmax;
       av_peakCa = av_peakCa + Camax;
       av_peakCafree = av_peakCafree + Cafreemax;
       w = w+1; p = p+1; s = s+1; q = q+1;
      
   else
  
       neff = neff-1; 
       av_tM = av_tM + tM;
       
   end
   
     
end
%end
%}
%% Shift and average

T12 = 30;

Is = zeros(1,N);

summ1=0;

if neff >0 
for m1=1:neff
Is = zeros(1,N);
I = Iset1(m1,:);

kin = round(t12(m1) - T12)/tstep;

if kin>=0
    for k=1:N-kin
        Is(k) = I(k+kin); 
    end
else
    kin = -kin;
    I(1:kin) = 0;
    for k=1+kin:N
        Is(k) = I(k-kin);
    end
end

Iset2(m1,:) = Is;

summ1=summ1+m1;

end

 if summ1>1
     Q = mean(Iset2);
 else
     Q = Is;
 end

else

Q = Is;
    
end

latencies = L;

%end
% Experimental results for QBs from Hardie
bump = load('bump_data0.txt'); 
QQ = bump(:,1);
num2 = length(QQ); % 65;
TT = 0:1:(num2-1);
T = 0:tstep:(num2-1);
QB = spline(TT,QQ,T);
Qmax = max(-QB);
T1 = 0;
T2 = 0;
result = zeros(1,num2); 

for k=1:length(T)
    if -QB(k) >= Qmax/2 && T1 == 0
        T1 = (k-1)*tstep;
    end
    if -QB(k) <= Qmax/2 && T1>0 && T2 == 0
        T2 = (k-1)*tstep;
    end
end

  TQ12 = (T1+T2)/2; 
  Tend = (length(T)-1)*tstep;
    Iav = Q;
        t1 = 0;
        t2 = 0; 
        [Imax,k50] = max(-Iav);
        for k=1:length(t)
            if -Iav(k) >= Imax/2 && t1 == 0
                t1 = (k-1)*tstep;
            end
            if -Iav(k) <= Imax/2 && k>k50 && t2 == 0
                t2 = (k-1)*tstep;
            end
        end
    t12 = (t1+t2)/2;
    istep = 1/tstep;
    kin = round(t12 - TQ12)/tstep; 
    kend = 1 + Tend/tstep;

    Is=zeros(1,kend);
 %       for k=1:kend
 %           Is(k) = Iav(k+kin); 
 %       end
if kin>=0
    for k=1:kend
        Is(k) = Iav(k+kin); 
    end
else
    kin = -kin;
    Iav(1:kin) = 0;
    for k=1+kin:kend
        Is(k) = Iav(k-kin);
    end
end

for i2 = 1:num2
    ir = (i2-1)*istep + 1;
    result(i2) = Is(ir);  % reduce the number of points to exp
end

%numel(T)
%numel(Q)

figure(10)
hold on
plot(TT,QQ,'or',T,Is,'k',TT,result,'x');
title(['Average QB (no.of bumps=',num2str(neff),') and experimental values (red circles)'])
xlabel('time [ms]')
ylabel('current [pA]')
hold off

av_Imax = mean(Imaxset);
sd_Imax = sqrt(cov(Imaxset));

av_L = mean(L);
sd_L = sqrt(cov(L));
av_t1 = mean(time1);
sd_t1 = sqrt(cov(time1));
av_t2 = mean(time2);
sd_t2 = sqrt(cov(time2));
av_t3 = mean(time3);
sd_t3 = sqrt(cov(time3));

av_tM = av_tM/runs;
av_peakPLC = av_peakPLC/neff;
av_PLCtot = av_PLCtot/neff;
av_peakDAG = av_peakDAG/neff;
av_peakNtrp = av_peakNtrp/neff;
av_peakCa = av_peakCa/neff;
%%av_peakCafree = av_peakCafree/neff;

format bank
format compact
disp('        tlat   t1=10 (sd=5.76)  t2=7.1 (2.6)  t3=12.85 (3.8)  Jpk=9.51 (3.8) ')
disp([av_L av_t1 av_t2 av_t3 av_Imax])
disp([sd_L sd_t1 sd_t2 sd_t3 sd_Imax])
disp(' ')
disp('        av_tM       av_peakPLC      av_PLCtot    av_peakDAG    av_peakNtrp    av_peakCa')
disp([av_tM av_peakPLC av_PLCtot av_peakDAG av_peakNtrp av_peakCa])
%av_tM av_peakPLC av_PLCtot av_peakDAG av_peakNtrp av_peakCa

%{
% Latency distribution ====================================================
iendL = length(latencies);
xL = 0:10:100;
boxL=zeros(1,11);
r=0;
for i=1:iendL
    if latencies(i)<5-r
        boxL(1)=boxL(1)+1;
        elseif latencies(i)<15-r
        boxL(2)=boxL(2)+1;
        elseif latencies(i)<25-r
        boxL(3)=boxL(3)+1;
        elseif latencies(i)<35-r
        boxL(4)=boxL(4)+1;
        elseif latencies(i)<45-r
        boxL(5)=boxL(5)+1;
        elseif latencies(i)<55-r
        boxL(6)=boxL(6)+1;
        elseif latencies(i)<65-r
        boxL(7)=boxL(7)+1;
        elseif latencies(i)<75-r
        boxL(8)=boxL(8)+1;
        elseif latencies(i)<85-r
        boxL(9)=boxL(9)+1;
        elseif latencies(i)<95-r
        boxL(10)=boxL(10)+1;
        elseif latencies(i)<105-r
        boxL(11)=boxL(11)+1;
    end
end

figure(20)
bar(xL,boxL,1,'y')
%}
% ===================================================================
% Latency distribution + 1ms (1ms estimated Rhodopsin activation time)
iendL = length(latencies);
xL1 = 0:10:100;
boxL1=zeros(1,11);
r=1;
for i=1:iendL
    if latencies(i)<5-r
        boxL1(1)=boxL1(1)+1;
        elseif latencies(i)<15-r
        boxL1(2)=boxL1(2)+1;
        elseif latencies(i)<25-r
        boxL1(3)=boxL1(3)+1;
        elseif latencies(i)<35-r
        boxL1(4)=boxL1(4)+1;
        elseif latencies(i)<45-r
        boxL1(5)=boxL1(5)+1;
        elseif latencies(i)<55-r
        boxL1(6)=boxL1(6)+1;
        elseif latencies(i)<65-r
        boxL1(7)=boxL1(7)+1;
        elseif latencies(i)<75-r
        boxL1(8)=boxL1(8)+1;
        elseif latencies(i)<85-r
        boxL1(9)=boxL1(9)+1;
        elseif latencies(i)<95-r
        boxL1(10)=boxL1(10)+1;
        elseif latencies(i)<105-r
        boxL1(11)=boxL1(11)+1;
    end
end

figure(21)
bar(xL1,boxL1,1,'g')
title('Distribution of latency times')
xlabel('time [ms]')
%ylabel('number')

% ===============================================================
% Imax distribution
iend = length(Imaxall);
x = 1:2:24;
box=zeros(1,12);
for i=1:iend
%     if Imaxall(i)<1
%        box(1)=box(1)+1;
    if Imaxall(i)<2
        box(1)=box(1)+1;
 %       elseif Imaxall(i)<3
 %       box(3)=box(3)+1;
        elseif Imaxall(i)<4
        box(2)=box(2)+1;
 %       elseif Imaxall(i)<5
 %       box(5)=box(5)+1;
        elseif Imaxall(i)<6
        box(3)=box(3)+1;
 %       elseif Imaxall(i)<7
 %       box(7)=box(7)+1;
        elseif Imaxall(i)<8
        box(4)=box(4)+1;
 %       elseif Imaxall(i)<9
 %       box(9)=box(9)+1;
        elseif Imaxall(i)<10
        box(5)=box(5)+1;
 %       elseif Imaxall(i)<11
 %       box(11)=box(11)+1;
        elseif Imaxall(i)<12
        box(6)=box(6)+1;
  %      elseif Imaxall(i)<13
  %      box(13)=box(13)+1;
        elseif Imaxall(i)<14
        box(7)=box(7)+1;
  %      elseif Imaxall(i)<15
  %      box(15)=box(15)+1;
        elseif Imaxall(i)<16
        box(8)=box(8)+1;
  %      elseif Imaxall(i)<17
  %      box(17)=box(17)+1;
        elseif Imaxall(i)<18
        box(9)=box(9)+1;
  %      elseif Imaxall(i)<19
  %      box(19)=box(19)+1;
        elseif Imaxall(i)<20
        box(10)=box(10)+1;
       elseif Imaxall(i)<22
        box(11)=box(11)+1;
       elseif Imaxall(i)<24
        box(12)=box(12)+1;
    end
end

figure(14)
hold on
bar(x,box,1,'y')
title('Maximum current distribution')
xlabel('current peak [pA]')
%ylabel('number')

%=========================================================================
endtime=clock;
%matlabpool close
time_taken=etime(endtime,starttime)

%% Plot <QB>

%compare_new(real(Q),t,tstep);




