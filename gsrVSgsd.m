clc, clear all, close all, warning off
clf;
run_time = 15000;
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
gL=0.1; gNa=1.5; gK=2; gsd=0.25; 
vL=-60; 
SK=0.25; v0K=-25; vK=-90; 
Ssd=0.09; v0sd=-40; vsd=50; 
SNa=0.25; v0Na=-25; vNa=50; 
vsr=-90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gsr=0.4;
T0=298;


tic
j=0;
i=0;
for gsr=0.2:0.01:0.7
    j=j+1;
for gsd=0:0.01:0.6
    
i=i+1;
T=291;
Phi=3^((T-T0)/10);
RO=1.3^((T-T0)/10);
a1=1./(1+exp(-SK*(-60-v0K)));
b1=1./(1+exp(-Ssd*(-60-v0sd))); 
c1=-(n./k)*RO*gsd*b1*(-60-vsd);
init1 = [-60 a1 b1 c1];
Er = 1e-5;
% Er = 1e-20;
options = odeset('RelTol',Er,'AbsTol',[Er Er Er Er]);
[T,V] = ode113(@System,[0 run_time],init1,options);
v = V(:,1);
aK = V(:,2);
asd = V(:,3);
asr = V(:,4);
p = 0.5;
sizi = length(v);
v1 = v(floor(sizi*p):end,1);
aK1 = aK(floor(sizi*p):end,1);
asd1 = asd(floor(sizi*p):end,1);
asr1 = asr(floor(sizi*p):end,1);
T = T(floor(sizi*p):end,1);
[~,loc] = findpeaks(v1);
spikeTimes1 = T(loc)        % get times when spikes occurred (ms)
spikeIntervals1 = spikeTimes1(2:length(spikeTimes1)) - spikeTimes1(1:length(spikeTimes1) - 1);
ISI1=max(spikeIntervals1);
peaks1=findpeaks(v1);
dpeak1=peaks1+60;
dpmin1 =abs( min(dpeak1))
if dpmin1>5
f1(i,j)=1000./ISI1
else
    f1(i,j)=0;
end

      %%%%%%%%%%%%%%%%%%%%%%%%%% rang kardan noghat
      
%       hold on
%       sz = 15;        %size of circles
%       c =  f1(i,j);    %color of circles
%       scatter(gsr,gsd,sz,c,'filled');

% T=285;
% Phi=3^((T-T0)/10);
% RO=1.3^((T-T0)/10);
% a2=1./(1+exp(-SK*(-60-v0K)));
% b2=1./(1+exp(-Ssd*(-60-v0sd))); 
% c2=-(n./k)*RO*gsd*b2*(-60-vsd);
% init2 = [-60 a2 b2 c2];
% 
% [T,V] = ode113(@System,[0 run_time],init2,options);
% v = V(:,1);
% aK = V(:,2);
% asd = V(:,3);
% asr = V(:,4);
% p = 0.8;
% sizi = length(v);
% v2 = v(floor(sizi*p):end,1);
% aK2 = aK(floor(sizi*p):end,1);
% asd2 = asd(floor(sizi*p):end,1);
% asr2 = asr(floor(sizi*p):end,1);
% T = T(floor(sizi*p):end,1);
% [~,loc] = findpeaks(v2);
% spikeTimes2 = T(loc)        % get times when spikes occurred (ms)
% spikeIntervals2 = spikeTimes2(2:length(spikeTimes2)) - spikeTimes2(1:length(spikeTimes2) - 1);
% ISI2=max(spikeIntervals2);
% peaks2=findpeaks(v2);
% dpeak2=peaks2+60;
% dpmin2 =abs( min(dpeak2))
% if dpmin2>5
% f2(i,j)=1000./ISI2
% else
%     f2(i,j)=0;
% end
% 
% T=279;
% Phi=3^((T-T0)/10);
% RO=1.3^((T-T0)/10);
% a3=1./(1+exp(-SK*(-60-v0K)));
% b3=1./(1+exp(-Ssd*(-60-v0sd))); 
% c3=-(n./k)*RO*gsd*b3*(-60-vsd);
% init3 = [-60 a3 b3 c3];
% 
% [T,V] = ode113(@System,[0 run_time],init3,options);
% v = V(:,1);
% aK = V(:,2);
% asd = V(:,3);
% asr = V(:,4);
% p = 0.5;
% sizi = length(v);
% v3 = v(floor(sizi*p):end,1);
% aK3 = aK(floor(sizi*p):end,1);
% asd3 = asd(floor(sizi*p):end,1);
% asr3 = asr(floor(sizi*p):end,1);
% T = T(floor(sizi*p):end,1);
% [~,loc] = findpeaks(v3);
% spikeTimes3 = T(loc)        % get times when spikes occurred (ms)
% spikeIntervals3 = spikeTimes3(2:length(spikeTimes3)) - spikeTimes3(1:length(spikeTimes3) - 1);
% ISI3=max(spikeIntervals3);
% peaks3=findpeaks(v3);
% dpeak3=peaks3+60;
% dpmin3 =abs( min(dpeak3))
% if dpmin3>5
% f3(i,j)=1000./ISI3
% else
%     f3(i,j)=0;
% end



end
end
toc
ofoghi = 0:.01:0.6;                                
amodi = 0.2:.01:0.7;
surf(ofoghi,amodi,f1);
contourf(distance,50);
