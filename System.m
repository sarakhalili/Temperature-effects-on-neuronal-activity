function dM = System(t,M)

global Cm
global Phi
global RO
global Isyn
global TK Tsr Tsd n k
global gL gNa gK gsd gsr vL SK v0K vK Ssd v0sd vsd SNa v0Na vNa vsr

dM = zeros(4,1);    
v = M(1);
aK = M(2);
asd = M(3);
asr = M(4);

aKf=1/(1+exp(-SK*(v-v0K)));
asdf=1/(1+exp(-Ssd*(v-v0sd)));
aNa=1/(1+exp(-SNa*(v-v0Na)));

Il=gL*(v-vL);
INa=RO*gNa*aNa*(v-vNa);
IK=RO*gK*aK*(v-vK);
Isd=RO*gsd*asd*(v-vsd);
Isr=RO*gsr*asr*(v-vsr);

dM(1) = (-Il-INa-IK-Isd-Isr-Isyn)/Cm;
dM(2) = (Phi*(aKf-aK))/TK;
dM(3) = (Phi*(asdf-asd))/Tsd;
dM(4) = (-Phi*(n*Isd+k*asr))/Tsr;
