clc, clear all, close all, warning off

run_time = 9000;
global Cm
global Phi
global RO
global Isyn
global TK Tsr Tsd n k
global gL gNa gK gsd gsr vL SK v0K vK Ssd v0sd vsd SNa v0Na vNa vsr
Cm=1;
Isyn=0;
TK=2; Tsd=10; Tsr=20; 
n=0.012; k=0.17;
gL=0.1; gNa=1.5; gK=2;
vL=-60; 
SK=0.25; v0K=-25; vK=-90; 
Ssd=0.09; v0sd=-40; vsd=50; 
SNa=0.25; v0Na=-25; vNa=50; 
vsr=-90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gsr=0;
gsd=0.25; 
T0=298;
T1=291;
Phi=3^((T1-T0)/10);
RO=1.3^((T1-T0)/10);


a=1./(1+exp(-SK*(-60-v0K)));
b=1./(1+exp(-Ssd*(-60-v0sd))); 
c=-(n./k)*RO*gsd*b*(-60-vsd);

init = [-60 a b c];

Er = 1e-5;
% Er = 1e-20;
options = odeset('RelTol',Er,'AbsTol',[Er Er Er Er]);
[T,V] = ode113(@System,[0 run_time],init,options);


v = V(:,1);
aK = V(:,2);
asd = V(:,3);
asr = V(:,4);

p = 0.8;
sizi = length(v);
v = v(floor(sizi*p):end,1);
aK = aK(floor(sizi*p):end,1);
asd = asd(floor(sizi*p):end,1);
asr = asr(floor(sizi*p):end,1);
T = T(floor(sizi*p):end,1);

[~,loc] = findpeaks(v);
spikeTimes3 = T(loc);        % get times when spikes occurred (ms)
spikeIntervals3 = spikeTimes3(2:length(spikeTimes3)) - spikeTimes3(1:length(spikeTimes3) - 1);
ISI3=min(spikeIntervals3);
peaks3=findpeaks(v);
dpeak3=peaks3+60;
dpmin3 =abs( min(dpeak3));
if dpmin3>10
f3=1000./ISI3
else
    f3=0
end

figure(1)

plot(T,v,'r','LineWidth',2)
print_text1=['gsd=',num2str(gsd,4),'  ,gsr=',num2str(gsr,4),'  ,T=',num2str(T1,4)];
title(print_text1)

ylabel('Membrain Voltage(mV)','FontName','Times New Roman','fontsize',17)
xlabel('Time(mS)','FontName','Times New Roman','fontsize',17)
T_min=min(T);
T_max=max(T);
axis([T_min,T_max,-80,40])

