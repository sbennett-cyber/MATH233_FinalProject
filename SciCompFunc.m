function [out1,out2,out3] = SciCompFunc(Start, Mid, End, D, k, h,n)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
StpInX=((Mid-Start)/h)+1;
    InnerStp=StpInX-2;
    func1=@(x) 2*sin((pi*x)/2)-sin(pi*x)+4*sin(2*pi*x);
    func2=@(x) 2*sin((pi*x)/2)-sin(pi*x)+4*sin(2*pi*x);
%     func1=@(x) 1/2*sin((pi*x));
%     func2=@(x) sin((pi*x));
%      True1=@(x,t) 1/2*sin((pi*x))*exp(-((pi^2)/2)*t); True2=@(x,t)
%      sin((pi*x))*exp(-((pi^2)/2)*t);
    r=(D*k)/(2*h^2);
    r1=(D*k)/(h^2);
    BC=[0,0,0];
    u0=MtrxYnet_Setup(Start, Mid, End, StpInX, func1, func2);
    Plot(Start, End, StpInX, u0, BC, func2(Mid));
    %% Setup Matrix of coefficents
    A=sparse(MtrxYnet(r,InnerStp, 3)); %Takes in coefficent, number of unknown points per channel, and number of channels
    A1=sparse(MtrxYnet(r1,InnerStp, 3)); %Takes in coefficent, number of unknown points per channel, and number of channels
    B=speye(length(A),length(A));  
    LHS=B-A;
    RHS=B+A;
    B1=speye(length(A1),length(A1));  
    LHS1=B1-A1;
    u1=u0;
    InterFace1=0;
    i=0;
    while i<n
      %  u1=biconjgrad(LHS1,u1,u1*0+1);
        u1=conjgrad(LHS1,u1,u1*0+1);
       % [u1,flag]=bicgstabl(LHS1,u1);
 %      u1=LHS\(RHS*u1);
       % InterFace=(-1/9)*(u(InnerStp-1)+u(2*InnerStp-1)+u(2*InnerStp+2))+(4/9)*(u(InnerStp)+u(2*InnerStp)+u(2*InnerStp+1)); %%Change to the flux at that point
        InterFace1=(-1/9)*(u1(InnerStp-1)+u1(2*InnerStp-1)+u1(2*InnerStp+2))+(4/9)*(u1(InnerStp)+u1(2*InnerStp)+u1(2*InnerStp+1)); %%Change to the flux at that point
        Plot(Start, End, StpInX, u1, BC, InterFace1);
        
        %Plot(Start, End, StpInX, u, BC, InterFace1);

        i=i+k;
    end
    NNx=StpInX-2;
    out1=[0;u1(1:NNx);InterFace1];
    out2=[0;u1(NNx+1:2*NNx);InterFace1];
    out3=[InterFace1;u1(2*NNx+1:3*NNx);0];
end

