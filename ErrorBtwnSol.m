function [out] = ErrorBtwnSol(errorspaceVector,dx)
%%% Shayna Bennett - 04/16/2020
%%% This function finds the order of a method by comparing the solutions at
%%% the final time when either the spatial step or the timestep are varied. 
%%% This follows page 257 of LeVeque's "Finite Difference Methods for
%%% Ordinary and Partial Differential Equations".
%%% Inputs: 
%%%     errorspaceMatrix: A matrix whose columns are solutions at different
%%%     spatial steps or timesteps
%%%     NNx: Number of unknowns in each row
%%%     NNy: Number of rows
%%% Outputs: 
%%%     out: A vector containing the order
L=length(errorspaceVector);
OrderSpace=zeros(L-1,1);
for a=1:(L-1)
    OrderSpace(a,1)=log2((errorspaceVector(a))/(errorspaceVector(a+1)));
end
out=OrderSpace;
end
