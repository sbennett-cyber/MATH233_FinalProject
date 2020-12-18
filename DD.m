%% Crank-Nicolson for heat equation on a line with choice of Dirichlet or Flux boundary conditions
%% AUTHOR: Shayna Bennett
%% Date: 4/1/2019
%% Version: 6.00 -- THE REVISED 4 VERSION 
% Description: This code will use Crank-Nicolson to solve for the solution
% to the heat equation on a line as described in Leveque-Ch 9 pg 182-183.
% Note***: This will need to be altered to get the correct
% vector of boundary conditions for the particular problem.

% Help from https://hnagib.com/portfolio/heat-equation/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Values to input or change for each problem are denoted with ** or ***

%% Test problem:
%%Taken from internet
% 4du/dt=du/dxx 0<=x<=2
% u(0,t)=u(2,t)=0    
% u(x,0)=2sin(pix/2)-sin(pix)+4sin(2pix)
% alpha=1/4; i=50; L=2; n=20; k=1
%% Inputs:
% alpha - coefficent on the second derivative 
% h - change in space
% L - length of line
% n - length of time
% k - change in time
% g0 - vector of left-side boundary conditions for all time
% gL - vector of right-side boundary conditions for all time
% u - vector of spatial boundary at initial time

function [out]=DD(r, L, h, g01, g02, gL1, gL2, u) 
%% Set Values:

y=linspace(0,L,single((L/h)+1))'; %Location of spatial points - use single to prevent doubles that effect the number of points if L is a multiple of pi

i=length(y); %Number of spatial gridpoints

u0=u(2:end-1);
%% Dirichlet-Dirichlet
%%Create Left and Right matrix

off(1:i-3)=r;
onL(1:i-2)=1+2*r;
onR(1:i-2)=1-2*r;

L=diag(onL,0)+diag(-off,-1)+diag(-off,1); % Create the centered spatial difference matrix
R=diag(onR,0)+diag(off,-1)+diag(off,1); % Create the centered spatial difference matrix

    b=R*u0+[r*(g01+g02);zeros(i-4,1);r*(gL1+gL2)];

    u0=L\b; %%Solve Ax=b problem 
    Z=[g01;u0;gL1];
out=Z;