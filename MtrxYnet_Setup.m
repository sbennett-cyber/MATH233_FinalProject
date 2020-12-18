function [out] = MtrxYnet_Setup(Start, Mid, End, StpInx, func1, func2)
%UNTITLED2 Summary of this function goes here
%  This works for zero Dirichlet BC

x=linspace(Start, Mid, StpInx)';
x2=linspace(Mid, End, StpInx)';
u1=func1(x);
u2=func1(x);
u3=func2(x2);
out=[u1(2:end-1);u2(2:end-1);u3(2:end-1)];


end

