clear all
close all
clc

dx=[0.05,0.025,0.0125,0.00625];
dt=0.9*([0.05,0.025,0.0125,0.00625].^2);

True1=@(x,t) 2*sin(pi*x/2)*exp(-pi^2*t/16)-sin(pi*x)*exp(-pi^2*t/4)+4*sin(2*pi*x)*exp(-pi^2*t);
True2=@(x,t) 2*sin(pi*x/2)*exp(-pi^2*t/16)-sin(pi*x)*exp(-pi^2*t/4)+4*sin(2*pi*x)*exp(-pi^2*t);

D=1/4;
Start=0;
Mid=1;
End=2;
n=2;

[R1,Q1,P1] = SciCompFunc(Start, Mid, End, D, dt(1), dx(1),n);
[R2,Q2,P2] = SciCompFunc(Start, Mid, End, D, dt(2), dx(2),n);
[R3,Q3,P3] = SciCompFunc(Start, Mid, End, D, dt(3), dx(3),n);
[R4,Q4,P4] = SciCompFunc(Start, Mid, End, D, dt(4), dx(4),n);


%%  Error Between Solutions in Space

RR2=R1-R2(1:2:end);
RR3=R2-R3(1:2:end);
RR4=R3-R4(1:2:end);
NormsR0=[normSB1D(RR2,dx(1),0),normSB1D(RR3,dx(2),0),normSB1D(RR4,dx(3),0)];
NormsR2=[normSB1D(RR2,dx(1),2),normSB1D(RR3,dx(2),2),normSB1D(RR4,dx(3),2)];

OrderLine1_Inf=ErrorBtwnSol(NormsR0,1)
OrderLine1_Two=ErrorBtwnSol(NormsR2,1)

PP2=P2(1:2:end);
PP3=P3(1:2:end);
PP4=P4(1:2:end);
NormsP0=[normSB1D(P1-PP2,dx(1),0),normSB1D(P2-PP3,dx(2),0),normSB1D(P3-PP4,dx(3),0)];
NormsP2=[normSB1D(P1-PP2,dx(1),2),normSB1D(P2-PP3,dx(2),2),normSB1D(P3-PP4,dx(3),2)];

OrderLine3_Inf=ErrorBtwnSol(NormsP0,1)
OrderLine3_Two=ErrorBtwnSol(NormsP2,1)

QQ2=Q2(1:2:end);
QQ3=Q3(1:2:end);
QQ4=Q4(1:2:end);
NormsQ0=[normSB1D(Q1-QQ2,dx(1),0),normSB1D(Q2-QQ3,dx(2),0),normSB1D(Q3-QQ4,dx(3),0)];
NormsQ2=[normSB1D(Q1-QQ2,dx(1),2),normSB1D(Q2-QQ3,dx(2),2),normSB1D(Q3-QQ4,dx(3),2)];

OrderLine2_Inf=ErrorBtwnSol(NormsQ0,1)
OrderLine2_Two=ErrorBtwnSol(NormsQ2,1)

