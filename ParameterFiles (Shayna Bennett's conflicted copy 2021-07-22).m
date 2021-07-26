%% Generate and Save Parameters

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16];
poly(1).Y = [0 16 16];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16];
poly(2).Y = [0 16 0];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[0 16];[0 16]];%Include one set of verticies for each row of Boundary
plot([0,0],[16 16])

D = [1,1];
ka = [1,1];
kd = [1,1];


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


save('SL001616kd1_Net.mat', 'net')
save('SL001616kd1_poly.mat', 'poly')

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 8 8 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [8 8 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[8 8];[8 8]];%Include one set of verticies for each row of Boundary
plot([0,16],[8 8])

D = [1,1];
ka = [0,0];
kd = [0,0];


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


save('StraightLinekd0V1_Net.mat', 'net')
save('StraightLinekd0V1_poly.mat', 'poly')

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 8.1415 8.1415 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [8.1415 8.1415 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8.1415;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[8.1415 8.1415];[8.1415 8.1415]];%Include one set of verticies for each row of Boundary
plot([0,16],[8.1415 8.1415])

D = [1,1];
ka = [1,1];
kd = [1,1];


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


save('StraightLinekd1V2_Net.mat', 'net')
save('StraightLinekd1V2_poly.mat', 'poly')
clear all
close all

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 8.1415 8.1415 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [8.1415 8.1415 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8.1415;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[8.1415 8.1415];[8.1415 8.1415]];%Include one set of verticies for each row of Boundary
plot([0,16],[8.1415 8.1415])

D = [4,4];
ka = [1,1];
kd = [1,1];


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


save('StraightLinekd1V2D4_Net.mat', 'net')
save('StraightLinekd1V2D4_poly.mat', 'poly')
clear all
close all


%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 8.1415 8.1415 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [8.1415 8.1415 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8.1415;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[8.1415 8.1415];[8.1415 8.1415]];%Include one set of verticies for each row of Boundary
plot([0,16],[8.1415 8.1415])

D = [1,1];
ka = [0,0];
kd = [0,0];


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


save('StraightLinekd0V2_Net.mat', 'net')
save('StraightLinekd0V2_poly.mat', 'poly')
clear all
close all

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 5.14 16*0.1415+5.14 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [5.14 16*0.1415+5.14 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 0.1415*8+5.14;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[5.14 16*0.1415+5.14];[5.14 16*0.1415+5.14]];%Include one set of verticies for each row of Boundary
plot([0,16],[5.14 16*0.1415+5.14])

D = [1,1];
ka = [1,1];
kd = [1,1];


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


save('StraightLinekd1_Net.mat', 'net')
save('StraightLinekd1_poly.mat', 'poly')
clear all
close all

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 16-0.1415 16 16 0];
poly(1).Y = [0.1415 16 16 0 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16-0.1415 0 0];
poly(2).Y = [0.1415 16 16 0];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 8.1415;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16-0.1415];[0,16-0.1415]]; %Include one set of verticies for each row of Boundary
vertexY = [[0.1415 16];[0.1415 16]];%Include one set of verticies for each row of Boundary
plot([0,16-0.1415],[0.1415 16])

D = [1,1];
ka = [1,1];
kd = [1,1];


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


save('StraightLinekd1V3_Net.mat', 'net')
save('StraightLinekd1V3_poly.mat', 'poly')
clear all
close all

%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 5.14 16*0.1415+5.14 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [5.14 16*0.1415+5.14 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 0.1415*8+5.14;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[5.14 16*0.1415+5.14];[5.14 16*0.1415+5.14]];%Include one set of verticies for each row of Boundary
plot([0,16],[5.14 16*0.1415+5.14])

D = [4,4];
ka = [1,1];
kd = [1,1];


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


save('StraightLinekd1D4_Net.mat', 'net')
save('StraightLinekd1D4_poly.mat', 'poly')
clear all
close all
%% Setup Single Line kd = ka = 0
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc
poly(1).X = [0 0 16 16];
poly(1).Y = [0 5.14 16*0.1415+5.14 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 16 16 0];
poly(2).Y = [5.14 16*0.1415+5.14 16 16];
poly(2).D = 1;
poly(2).r = 1;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [1,2;2,1];
LineID = [1,1];
SharedVertexX = [];
SharedVertexY = [];
SharedWith = [1,2];
vertexX = [[0,16];[0,16]]; %Include one set of verticies for each row of Boundary
vertexY = [[5.14 16*0.1415+5.14];[5.14 16*0.1415+5.14]];%Include one set of verticies for each row of Boundary
plot([0,16],[5.14 16*0.1415+5.14])

D = [1,1];
ka = [0,0];
kd = [0,0];


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
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 0.1415*8+5.14;


save('StraightLinekd0_Net.mat', 'net')
save('StraightLinekd0_poly.mat', 'poly')
clear all
close all
%% Setup for Kinked Line Case 2
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc

poly(1).X = [0 10 16 16 0];
poly(1).Y = [0 10 0 16 16];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 10 16];
poly(2).Y = [0 10 0];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 10; %The Gaussian should be centered on the network
net.CenterGaussianY = 10;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[2,1];[2,1];[1,2]];
LineID = [1,2,1,2];
SharedVertexX = [10];
SharedVertexY = [10];
SharedWith = [1,2];
vertexX = [[0,10];[10,16];[0,10];[10,16]];
vertexY = [[0,10];[10,0];[0,10];[10,0]];
plot([0,10],[0,10])
hold on
plot([10,16],[10,0])

D = [1,1,1,1];
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

save('KinkedLine_Net.mat', 'net')
save('KinkedLine_poly.mat', 'poly')
clear all
close all

%% Setup for Kinked Line Case 2
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc

poly(1).X = [0 10.1415 16 16 0];
poly(1).Y = [4 10.1415 5 16 16];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 10.1415 16 16 0];
poly(2).Y = [4 10.1415 5 0 0];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 10.1415; %The Gaussian should be centered on the network
net.CenterGaussianY = 10.1415;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[2,1];[2,1];[1,2]];
LineID = [1,2,1,2];
SharedVertexX = [10.1415];
SharedVertexY = [10.1415];
SharedWith = [1,2];
vertexX = [[0,10.1415];[10.1415,16];[0,10.1415];[10.1415,16]];
vertexY = [[4,10.1415];[10.1415,5];[4,10.1415];[10.1415,5]];
plot([0,10.1415],[0,10.1415])
hold on
plot([10.1415,16],[10.1415,0])

D = [4,4,4,4];
% ka = [1,1,1,1];
% kd = [1,1,1,1];
ka = [0,0,0,0];
kd = [0,0,0,0];


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

save('KinkedLinekd0V2_Net.mat', 'net')
save('KinkedLinekd0V2_poly.mat', 'poly')
clear all
close all

%% Setup for Kinked Line Case 2
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc

poly(1).X = [0 10 16 16 0];
poly(1).Y = [0.1415 10.1415 0 16 16];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 10 16 0];
poly(2).Y = [0.1415 10.1415 0 0];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 10; %The Gaussian should be centered on the network
net.CenterGaussianY = 7;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[2,1];[2,1];[1,2]];
LineID = [1,2,1,2];
SharedVertexX = [10];
SharedVertexY = [10.1415];
SharedWith = [1,2];
vertexX = [[0,10];[10,16];[0,10];[10,16]];
vertexY = [[0.1415,10.1415];[10.1415,0];[0.1415,10.1415];[10.1415,0]];
plot([0,10],[0.1415,10.1415])
hold on
plot([10,16],[10.1415,0])

D = [1,1,1,1];
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

save('KinkedLinekd1V6_Net.mat', 'net')
save('KinkedLinekd1V6_poly.mat', 'poly')
clear all
close all

%% Setup for Kinked Line Case 2
% Setup the verteces as a structure [coordX, coordY]
clear all
close all
clc

poly(1).X = [0 0 10 16 16];
poly(1).Y = [0 6.1415 10.1415 6 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [0 10 16 16 0];
poly(2).Y = [6.1415 10.1415 6 16 16];
poly(2).D = 1;
poly(2).r = 1;
net.CenterGaussianX = 10; %The Gaussian should be centered on the network
net.CenterGaussianY = 10.1415;


% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[2,1];[2,1];[1,2]];
LineID = [1,2,1,2];
SharedVertexX = [10];
SharedVertexY = [10.1415];
SharedWith = [1,2];
vertexX = [[0,10];[10,16];[0,10];[10,16]];
vertexY = [[6.1415,10.1415];[10.1415,6];[6.1415,10.1415];[10.1415,6]];
plot([0,10],[6.1415,10.1415])
hold on
plot([10,16],[10.1415,6])

D = [1,1,1,1];
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

save('KinkedLinekd1V7_Net.mat', 'net')
save('KinkedLinekd1V7_poly.mat', 'poly')
clear all
close all

%% Setup for Ynet Case 2
% Setup the verteces as a structure [coordX, coordY]

poly(1).X = [0 0 8 16];
poly(1).Y = [0 8 8 0];
poly(1).D = 1;
poly(1).r = 1;

poly(2).X = [8 16 16];
poly(2).Y = [8 0 16];
poly(2).D = 1;
poly(2).r = 1;

poly(3).X = [0 8 16 0];
poly(3).Y = [8 8 16 16];
poly(3).D = 1;
poly(3).r = 1;
net.CenterGaussianX = 8; %The Gaussian should be centered on the network
net.CenterGaussianY = 6;

% Setup the network as a structure [[Domain of origin, Cross], LineID, coordX, coordY, ka, kd, D]
Boundary = [[1,2];[1,3];[2,1];[2,3];[3,1];[3,2]]; 
LineID = [2,1,2,3,1,3];
SharedVertexX = [10];
SharedVertexY = [16];
SharedWith = [1,2,3];
vertexX = [[8,16];[0,8];[8,16];[8,16];[0,8];[8,16]];
vertexY = [[8,0];[8,8];[8,0];[8,16];[8,8];[8,16]];
plot([8,16],[8,0])
hold on
plot([0,8],[8,8])
hold on
plot([8,16],[8,16])

D = [1,1,1,1,1,1];
ka = [1,1,1,1,1,1];
kd = [1,1,1,1,1,1];


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

save('Ynet_Net.mat', 'net')
save('Ynet_Poly.mat', 'poly')