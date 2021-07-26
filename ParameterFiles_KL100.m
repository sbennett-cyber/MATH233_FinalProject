%% Setup for Kinked Line Case 2
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc

poly(1).X = [0 100 100 50.1415 0];
poly(1).Y = [0 0 50 60.1415 50];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 50.1415 100 100 0];
poly(2).Y = [50 60.1415 50 100 100];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 50; %The Gaussian should be centered on the network
net.CenterGaussianY = 50;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[2,1];[2,1];[1,2]];
LineID = [1,2,1,2];
SharedVertexX = [50.1415];
SharedVertexY = [60.1415];
SharedWith = [1,2];
% vertexX = [[0,50.1415];[50.1415,100];[0,50.1415];[50.1415,100]];
% vertexY = [[50,60.1415];[60.1415,50];[50,60.1415];[60.1415,50]];
vertexX = [[0,50.1415];[50.1415,100];[0,50.1415];[50.1415,100]];
vertexY = [[50,60.1415];[60.1415,50];[50,60.1415];[60.1415,50]];
plot([0,50.1415],[50,60.1415])
hold on
plot([50.1415,100],[60.1415,50])

D = [8,8,8,8];
ka = [1,1,1,1];
kd = [1,1,1,1];


for i = 1:2
    net.Boundary(:,i) = Boundary(:,i); 
    net.vertexX(:,i) = vertexX(:,i);
    net.vertexY(:,i) = vertexY(:,i);
end    
net.LineID = LineID;
net.D = D;
net.ka = ka;
net.kd = kd;
net.SharedVertexX = SharedVertexX;
net.SharedVertexY = SharedVertexY;
net.SharedWith = SharedWith;

save('KL00100100kd1D8V1_Net.mat', 'net')
save('KL00100100kd1D8V1_poly.mat', 'poly')
clear all
close all