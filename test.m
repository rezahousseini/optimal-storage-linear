clear all;
close all;

filename="optimalStorage";

numN=100;
T=1/4;
n=0:numN-1;
t=n*T;

S.Qmax=[7,Inf];
S.Qmin=[0,Inf];
S.q0=[7,Inf];
S.C=[1,10];
S.D=[1,10];
S.nul=[1,1];
S.nuc=[1,1];
S.nud=[1,1];
S.DeltaCmax=[0.1,10];
S.DeltaDmax=[0.2,10];

numS=length(S.Qmax);

u=zeros(1,numN);
for k= 2:numN
	u(k)=0.9*u(k-1)+normrnd(0,0.1,1,1);
endfor

g=exp(0.16+0.1*cos(2*pi.*n*T/24-5*pi/4)+u+normrnd(0,0.1,1,numN));
p=ones(1,numN);%exp(0.1+0.4*cos(2*pi.*n*T/24-3*pi/2)+u+normrnd(0,0.1,1,numN));
r=exp(0.18+0.3*cos(2*pi.*n*T/24-pi)+u+normrnd(0,0.1,1,numN));

% Set of time steps
Set.I=1:numN;
% Set of storage variables with Qmax < inf
Set.F=[1];
% Set of storage variables with Qmax = Inf
Set.E=[2];

numSfin=length(Set.F);

Param.N=numN;
Param.g=g;
Param.r=r;
Param.pg=zeros(1,numN);
Param.pr=zeros(1,numN);
Param.pc=[zeros(numS-1,numN);p];
Param.pd=[zeros(numS-1,numN);p];
Param.T=T;
Param.Qmax=S.Qmax;
Param.Qmin=S.Qmin;
Param.q0=S.q0;
Param.C=S.C;
Param.D=S.D;
Param.nul=S.nul;
Param.nuc=S.nuc;
Param.nud=S.nud;
Param.DeltaCmax=S.DeltaCmax;
Param.DeltaDmax=S.DeltaDmax;

%mkoctfile writeOptimizationModelData.cc -lglpk

writeOptimizationModelData(filename,Set,Param);

%mkoctfile optimalNStorage.cc -lglpk

tic;
[q,uc,ud,cost]=optimalNStorage(filename,Set,Param);
toc;

figure(1)
plot(t,r,t,g,t,g+ud(1,:)-uc(1,:),t,ud(2,:)-uc(2,:))
legend("Nachfrage","Produktion ohne Speicher","Produktion mit Speicher","Netz")
