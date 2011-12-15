clear all;
close all;

filename="optimalStorage";

N=100;
T=1;
n=0:N-1;
t=n*T;

A=2;

S.Qmax=[1,Inf];
S.Qmin=[0,0];
S.q0=[1,Inf];
S.C=[0.3,Inf];
S.D=[0.3,Inf];
S.nul=[1,1];
S.nuc=[1,1];
S.nud=[1,1];

u=zeros(1,N);
for k= 2:N
	u(k)=0.9*u(k-1)+normrnd(0,0.1,1,1);
endfor

g=exp(0.16+0.1*cos(2*pi.*n*T/24-5*pi/4)+u+normrnd(0,0.1,1,N));
p=exp(0.1+0.4*cos(2*pi.*n*T/24-3*pi/2)+u+normrnd(0,0.1,1,N));
r=exp(0.18+0.4*cos(2*pi.*n*T/24-pi)+u+normrnd(0,0.1,1,N));

% Set of time steps
Set.I=1:N;
% Set of storage variables with no inf
Set.A=1;
% Set of storage variables with C,D and Qmax = Inf
Set.E=2;
% Set of storage variables with Qmax = Inf
Set.F=[];

Param.N=N;
Param.g=g;
Param.r=r;
Param.pg=zeros(1,N);
Param.pr=zeros(1,N);
Param.pc=[zeros(1,N);p];
Param.pd=zeros(A,N);
Param.T=T;
Param.Qmax=S.Qmax;
Param.Qmin=S.Qmin;
Param.q0=S.q0;
Param.C=S.C;
Param.D=S.D;
Param.nul=S.nul;
Param.nuc=S.nuc;
Param.nud=S.nud;

%mkoctfile writeOptimizationModelData.cc -lglpk

writeOptimizationModelData(filename,Set,Param);

%mkoctfile optimal1Storage.cc -lglpk

tic;
[ucS,udS,ucG,udG]=optimal1Storage(filename,N);
toc;

%q=zeros(1,N);
%q(1)=q0(1);
%for k=2:N
%q(k)=nul(1)*q(k-1)+T*(nuc(1)*ucS(k))
%endfor

figure(1)
plot(t,r,t,g,t,g+udS-ucS)
legend("Nachfrage","Produktion ohne Speicher","Produktion mit Speicher")

%figure(2)
%plot(t,ucS-udS+ucG-udG+r-g)
