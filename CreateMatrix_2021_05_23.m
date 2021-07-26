% 
% % %% Test the use of polygons to define subdomains in my code
% % % 
% % % %Create the shape of the domain
% % % 
% % % pgon1 = polyshape([0 0 8 16],[0 8 8 0]);
% % % plot(pgon1)
% % % hold on
% % % pgon2 = polyshape([0 8 16 0],[8 8 16 16]);
% % % plot(pgon2)
% % % hold on
% % % pgon3 = polyshape([8 16 16],[8 0 16]);
% % % plot(pgon3)
% % % hold on
% % % 
% % % %Identify grid points within and outside of the subdomains (Identity function)
% % % X=linspace(0,16,11);
% % % [Xm,Ym]=meshgrid(X,X);
% % % xm=reshape(Xm,11^2,1);
% % % ym=reshape(Ym,11^2,1);
% % % 
% % % plgn1X=[0 0 8 16];
% % % plgn1Y= [0 8 8 0];
% % % [in1,out1]=inpolygon(xm,ym,plgn1X,plgn1Y);
% % % 
% % % plot(xm(in1),ym(in1),'r*')
% % % 
% % % normalX = [6,6];
% % % normalY = [9,7];
% % % 
% % % [xi,yi] = inpolygon(normalX,normalY,plgn1X,plgn1Y);
% % % 
% % % figure
% % % plot(pgon1)
% % % hold on
% % % plot(normalX,normalY)
% % % 
% % %% Tests for actual Code
% % 
% %% Updated to allow the EBM to work with all angles for KinkedLineV1, KinkedLineV2, and YnetV1
% 
clear all
close all
clc
% 
% Solve

load('Ynet_Net.mat')
load('Ynet_Poly.mat')
% load('KinkedLinekd0V2_Net.mat')
% load('KinkedLinekd0V2_poly.mat')
% load('StraightLinekd0V1_Net.mat')
% load('StraightLinekd0V1_poly.mat')
% load('SL00100100kd1D1V1_Net.mat')
% load('SL00100100kd1D1V1_poly.mat')
% load('SL00100100kd1_Net.mat')
% load('SL00100100kd1_poly.mat')
% load('KL00100100A2kd1D8_Net.mat')
% load('KL00100100A2kd1D8_poly.mat')

% Setup discretization points
% PntsInX = 251;
% PntsInY = 251;
% PntsLine = [501, 0, 0]; %This is a vector when the network has more than one channel
% PntsInX = 201;
% PntsInY = 201;
% PntsLine = [401, 401, 0]; %This is a vector when the network has more than one channel

PntsInX = 81;
PntsInY = 81;
PntsLine = [161, 161, 161]; %This is a vector when the network has more than one channel

alpha = 1;
beta = 1;
%CenterGaussianY = 8.1415;
L = 16;
dt=0.2;
n = 20;
Time = (n/dt);
C = NetworkLine(poly,net, [0,L],[0,L],[PntsInX,PntsInY], PntsLine, alpha, beta); %Initialize the class with polygons (vertec of each subdomain & the diffusion and growth terms) and net from .mat files, min and max in X, min and max in Y, [Pnts in X,Pnts in y,dt]
Poly = C.CreatePolygon(poly, "Plot"); %Initialize subdomains and returns structure [coordX, coordY, polygon, PntsInX, PntsInY]. The last two are the coordinates of the gridpoints that fall in the domain

NetworkMatrix = C.CreateMatrixNetwork(dt/2,"ZeroDirichlet");
FieldMatrix = C.CreateMatrixFields(Poly,dt,"ZeroDirichlet", "NoPlot");

[Line, FieldIC] = C.SetIC(Poly,NetworkMatrix.RhoID_X ,NetworkMatrix.RhoID_Y);

C.PlotSolutions(Poly, FieldIC)

i = 0
while i < Time

    Line1 = C.SolveNetwork(dt/2, Poly, NetworkMatrix ,Line, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
    FieldIC  = C.SolveField (i, dt, Poly, NetworkMatrix ,FieldMatrix ,Line1, FieldIC, "Solution"); %Use Line if studying the Laplacian and Line1 for the full problem
    
    Line = C.SolveNetwork(dt/2, Poly, NetworkMatrix ,Line1, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
    %figure
    C.PlotSolutions(Poly, FieldIC)
    %figure
    %plot (NetworkMatrix.RhoID_X, Line.Line)
    %hold on
    drawnow
    i = i + 1;
end


%% Convergence Study = Laplacian
clear all
close all
clc

Conv = ImportantFunctions();
load('StraightLinekd1V1_Net.mat')
load('StraightLinekd1V1_poly.mat')
% load('KinkedLineV2_Net.mat')
% load('KinkedLineV2_poly.mat')
alpha = 1;
beta = 1;
PntsInX  = [21; 41; 81; 161; 321];
PntsLine = [41,0,0; 81,0,0; 161,0,0; 321,0,0; 641,0,0];
% drho = 16./(PntsLine(:,1)-1);
% dt=[1,1,1,1,1];
% dtStep = [1,2,4,8,16];

% PntsInX  = [21; 41; 81; 161; 321];
% PntsLine = [41,41,0; 81,81,0;161,161,0; 321,321,0; 641,641,0];
drho = 16./(PntsLine(:,1)-1);
dt=[1; 1; 1; 1; 1];
dx = 16./(PntsInX(:,1)-1);
dtStep = [1,1,1,1,1];

for i = 1:length(PntsInX)
    C = NetworkLine(poly,net, [0,16],[0,16],[PntsInX(i),PntsInX(i)], PntsLine(i,:), alpha, beta); %Initialize the class with polygons (vertec of each subdomain & the diffusion and growth terms) and net from .mat files, min and max in X, min and max in Y, [Pnts in X,Pnts in y,dt]
    Poly = C.CreatePolygon(poly, "NoPlot"); %Initialize subdomains and returns structure [coordX, coordY, polygon, PntsInX, PntsInY]. The last two are the coordinates of the gridpoints that fall in the domain

    NetworkMatrix = C.CreateMatrixNetwork(dt(i),"ZeroDirichlet");
    FieldMatrix = C.CreateMatrixFields(Poly,dt(i),"ZeroDirichlet", "NoPlot");

    [Line, FieldIC] = C.SetIC(Poly,NetworkMatrix.RhoID_X ,NetworkMatrix.RhoID_Y);
    
    figure
    plot (NetworkMatrix.RhoID_X, Line.Line)
    hold on
    drawnow
    FieldNB = sparse((PntsInX(i)-2)^2-1);
    Field = sparse(PntsInX(i),PntsInX(i));
    for j = 1:dtStep(i)
     %   Line1 = C.SolveNetwork(dt(j)/2, Poly, NetworkMatrix ,Line, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
        FieldIC = C.SolveField (j, dt(j), Poly, NetworkMatrix ,FieldMatrix ,Line, FieldIC, "LaplaceError");
    
    %    Line = C.SolveNetwork(dt(j)/2, Poly, NetworkMatrix ,Line1, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    C.PlotSolutions(Poly, FieldIC);
%     figure
%     plot (NetworkMatrix.RhoID_X, Line)
    hold on
    drawnow
    
    end
%     for k = 1:2
%         [in,bd] = inpolygon(Poly(k).InXNB,Poly(k).InYNB,Poly(k).X,Poly(k).Y);
%         PntsOnNet = find(bd==1);
%         FieldIC(k).Field(PntsOnNet) = 0;
%     end
    FieldNB(Poly(2).LocateNB) = FieldIC(2).Field;
    FieldNB(Poly(1).LocateNB) = FieldIC(1).Field;
    FieldNB2 = reshape(FieldNB,PntsInX(i)-2,PntsInX(i)-2)';
    Field(2:end-1,2:end-1) = FieldNB2;
    Results(i).Line(:,:) = Field;
end

[TB_L] = Conv.ErrorWTrueSol2D(Results(1).Line, Results(2).Line, Results(3).Line, Results(4).Line, Results(5).Line,dx,dt)
%    
% %I still need to finish the function to create network matrix for three
% %channels, finish the function to solve in time for concentration on two
% %and three channels, make poly a property of the class, write function to
% %interpolate the concentration of rho onto the C_Gamma points. 
% 
%% Convergence Study = Both
clear all
close all
clc

Conv = ImportantFunctions();

load('KinkedLinekd1V3_Net.mat')
load('KinkedLinekd1V3_poly.mat')
% load('StraightLinekd1V2_Net.mat')
% load('StraightLinekd1V2_poly.mat')
% load('SL00100100kd1_Net.mat')
% load('SL00100100kd1_poly.mat')

alpha = 1;
beta = 1;
n = 20;
L = 16;
PntsInX  = [21; 41; 81; 161; 321];
PntsLine = [41,41,0; 81,81,0;161,161,0; 321,321,0; 641,641,0];
% PntsInX  = [101; 201; 401; 801; 1601];
% PntsLine = [201,201,0; 401,401,0;801,801,0; 1601,1601,0; 3201,3201,0];
drho = L./(PntsLine(:,1)-1);
dt = L./(PntsInX(:,1)-1);
dx = L./(PntsInX(:,1)-1);


for i = 1:length(PntsInX)
    tic
    C = NetworkLine(poly,net, [0,L],[0,L],[PntsInX(i),PntsInX(i)], PntsLine(i,:), alpha, beta); %Initialize the class with polygons (vertec of each subdomain & the diffusion and growth terms) and net from .mat files, min and max in X, min and max in Y, [Pnts in X,Pnts in y,dt]
    Poly = C.CreatePolygon(poly, "Plot"); %Initialize subdomains and returns structure [coordX, coordY, polygon, PntsInX, PntsInY]. The last two are the coordinates of the gridpoints that fall in the domain

    NetworkMatrix = C.CreateMatrixNetwork(dt(i)/2,"ZeroDirichlet");
    FieldMatrix = C.CreateMatrixFields(Poly,dt(i),"ZeroDirichlet", "NoPlot");

    [Line, FieldIC] = C.SetIC(Poly,NetworkMatrix.RhoID_X ,NetworkMatrix.RhoID_Y);
    
    figure
    plot (NetworkMatrix.RhoID_X, Line.Line)
    hold on
    drawnow
    FieldNB = sparse((PntsInX(i)-2)^2-1);
    Field = sparse(PntsInX(i),PntsInX(i));
    
    figure
    j = 0;
    while j < n/dt(i)
        Line1 = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,Line, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
        FieldIC = C.SolveField (j, dt(i), Poly, NetworkMatrix ,FieldMatrix ,Line1, FieldIC, "Solution");
    
        Line = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,Line1, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
        
        j = j + 1;
%    C.PlotSolutions(Poly, FieldIC);
%    plot (NetworkMatrix.RhoID_X, Line.Line)
%    hold on
%    drawnow
    
    end
% %     for k = 1:2
% %         [in,bd] = inpolygon(Poly(k).InXNB,Poly(k).InYNB,Poly(k).X,Poly(k).Y);
% %         PntsOnNet = find(bd==1);
% %         FieldIC(k).Field(PntsOnNet) = 0;
% %     end
    
    FieldNB(Poly(1).LocateNB) = FieldIC(1).Field;
    FieldNB(Poly(2).LocateNB) = FieldIC(2).Field;
    FieldNB2 = reshape(FieldNB,PntsInX(i)-2,PntsInX(i)-2)';
    Field(2:end-1,2:end-1) = FieldNB2;
    Results(i).Field(:,:) = Field;
    Results(i).Line(:,:) = [0;Line.Line;0];
    Results(i).dx = dx;
    Results(i).drho = drho;
    Results(i).dt = dt;
    Results(i).FinalTime = n;
    toc
    i
end
%save('OutputStraightLinekd1V1_Both.mat','Results')
[TField_B] = Conv.ErrorBtwnSol2D(Results(1).Field, Results(2).Field, Results(3).Field, Results(4).Field, Results(5).Field,dx,dt,'Space')
% table2latex(TField_B,'StraightLinekd1V1_Field_Both')

[TLine_B]  = Conv.ErrorBtwnSol(Results(1).Line, Results(2).Line, Results(3).Line, Results(4).Line, Results(5).Line,drho,dt,'Space')
% table2latex(TLine_B,'StraightLinekd1V1_Line_Both')

%% Convergence Study = Space
clear all
close all
%clc

Conv = ImportantFunctions();

load('StraightLinekd1V1_Net.mat')
load('StraightLinekd1V1_poly.mat')
% load('SL00100100kd1_Net.mat')
% load('SL00100100kd1_poly.mat')

% load('SL001616kd1_Net.mat')
% load('SL001616kd1_poly.mat')

alpha = 1;
beta = 1;
n = 1.6;
L = 16;
% PntsInX  = [21; 41; 81; 161; 321];
% PntsLine = [41,0,0; 81,0,0;161,0,0; 321,0,0; 641,0,0];
PntsInX  = [21;41; 81; 161; 321];
PntsLine = [41,0,0;81,0,0;161,0,0; 321,0,0; 641,0,0];
drho = L./(PntsLine(:,1)-1);
dx = L./(PntsInX(:,1)-1);

dt = [0.0125;0.0125;0.0125;0.0125;0.0125];

for i = 1:length(PntsInX)
    tic
    C = NetworkLine(poly,net, [0,L],[0,L],[PntsInX(i),PntsInX(i)], PntsLine(i,:), alpha, beta); %Initialize the class with polygons (vertec of each subdomain & the diffusion and growth terms) and net from .mat files, min and max in X, min and max in Y, [Pnts in X,Pnts in y,dt]
    Poly = C.CreatePolygon(poly, "NoPlot"); %Initialize subdomains and returns structure [coordX, coordY, polygon, PntsInX, PntsInY]. The last two are the coordinates of the gridpoints that fall in the domain

    NetworkMatrix = C.CreateMatrixNetwork(dt(i)/2,"ZeroDirichlet");
    FieldMatrix = C.CreateMatrixFields(Poly,dt(i),"ZeroDirichlet", "NoPlot");

    [Line, FieldIC] = C.SetIC(Poly,NetworkMatrix.RhoID_X ,NetworkMatrix.RhoID_Y);

    figure
    plot (NetworkMatrix.RhoID_X, Line.Line)
    hold on
    drawnow
    FieldNB = sparse((PntsInX(i)-2)^2-1);
    Field = sparse(PntsInX(i),PntsInX(i));
    
    j=0;
    while j < n/dt(i)
        Line1 = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,Line, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
        FieldIC = C.SolveField (j, dt(i), Poly, NetworkMatrix ,FieldMatrix ,Line1, FieldIC, "Solution");
    
        Line = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,Line1, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
        j = j + 1;
    C.PlotSolutions(Poly, FieldIC);
%     figure
%     plot (NetworkMatrix.RhoID_X, Line)
    hold on
    drawnow
    
    end
%     for k = 1:2
%         [in,bd] = inpolygon(Poly(k).InXNB,Poly(k).InYNB,Poly(k).X,Poly(k).Y);
%         PntsOnNet = find(bd==1);
%         FieldIC(k).Field(PntsOnNet) = 0;
%     end
    toc
    FieldNB(Poly(1).LocateNB) = FieldIC(1).Field;
    FieldNB(Poly(2).LocateNB) = FieldIC(2).Field;
    FieldNB2 = reshape(FieldNB,PntsInX(i)-2,PntsInX(i)-2)';
    Field(2:end-1,2:end-1) = FieldNB2;
    Results(i).Field(:,:) = Field;
    Results(i).Line(:,:) = [0;Line.Line;0];
    Results(i).dx = dx;
    Results(i).drho = drho;
    Results(i).dt = dt;
    Results(i).FinalTime = n;
    toc
    i
end
%save('OutputSL00100100kd1_Space.mat','Results')
[TField_S] = Conv.ErrorBtwnSol2D(Results(1).Field, Results(2).Field, Results(3).Field, Results(4).Field, Results(5).Field,dx,dt,'Space')
%table2latex(TField_S,'SL00100100kd1_Field_Space')

[TLine_S]  = Conv.ErrorBtwnSol(Results(1).Line, Results(2).Line, Results(3).Line, Results(4).Line, Results(5).Line,drho,dt,'Space')
%table2latex(TLine_S,'SL00100100kd1_Line_Space')

%% Convergence Study = Time
clear all
close all
%clc

Conv = ImportantFunctions();
% load('StraightLinekd1_Net.mat')
% load('StraightLinekd1_poly.mat')
% load('SL001616kd1_Net.mat')
% load('SL001616kd1_poly.mat')
load('SL00100100kd1_Net.mat')
load('SL00100100kd1_poly.mat')

alpha = 1;
beta = 1;
n = 10;
L = 100;
PntsInX  = [641,641,641,641,641];
PntsLine = [1281,0,0; 1281,0,0;1281,0,0; 1281,0,0; 1281,0,0];%[41,0,0; 81,0,0;161,0,0; 321,0,0; 641,0,0];
drho = L./(PntsLine(:,1)-1);
dx = (L./(PntsInX-1))';


dt = [0.8;0.4;0.2;0.1;0.05];
dtStep = n./dt;
for i = 1:length(PntsInX)
    tic
    C = NetworkLine(poly,net, [0,L],[0,L],[PntsInX(i),PntsInX(i)], PntsLine(i,:), alpha, beta); %Initialize the class with polygons (vertec of each subdomain & the diffusion and growth terms) and net from .mat files, min and max in X, min and max in Y, [Pnts in X,Pnts in y,dt]
    Poly = C.CreatePolygon(poly, "NoPlot"); %Initialize subdomains and returns structure [coordX, coordY, polygon, PntsInX, PntsInY]. The last two are the coordinates of the gridpoints that fall in the domain

    NetworkMatrix = C.CreateMatrixNetwork(dt(i)/2,"ZeroDirichlet");
    FieldMatrix = C.CreateMatrixFields(Poly,dt(i),"ZeroDirichlet", "NoPlot");

    [LineIC, FieldIC] = C.SetIC(Poly,NetworkMatrix.RhoID_X ,NetworkMatrix.RhoID_Y);

    figure
    plot (NetworkMatrix.RhoID_X, LineIC.Line)
    hold on
    drawnow
    FieldNB = sparse((PntsInX(i)-2)^2-1);
    Field = sparse(PntsInX(i),PntsInX(i));
    
    j=0;
    
    while j < n/dt(i)
        Line1 = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,LineIC, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
    
        FieldIC = C.SolveField (j, dt(i), Poly, NetworkMatrix ,FieldMatrix ,Line1, FieldIC, "Solution");
    
        Line = C.SolveNetwork(dt(i)/2, Poly, NetworkMatrix ,Line1, FieldIC, "ZeroDirichlet", "NoPlot"); %The dt must be the same as used in the matrix
        j = j+1;
    end
    for k = 1:2
        [in,bd] = inpolygon(Poly(k).InXNB,Poly(k).InYNB,Poly(k).X,Poly(k).Y);
        PntsOnNet = find(bd==1);
        FieldIC(k).Field(PntsOnNet) = 0;
    end
    toc
    FieldNB(Poly(1).LocateNB) = FieldIC(1).Field;
    FieldNB(Poly(2).LocateNB) = FieldIC(2).Field;
    FieldNB2 = reshape(FieldNB,PntsInX(i)-2,PntsInX(i)-2)';
    Field(2:end-1,2:end-1) = FieldNB2;
    Results(i).Field(:,:) = Field;
    Results(i).Line(:,:) = [0;Line.Line;0];
    Results(i).dx = dx;
    Results(i).drho = drho;
    Results(i).dt = dt;
    Results(i).FinalTime = n;
    toc
    i
end
save('OutputSL00100100kd1_Time.mat','Results')
[TField_T] = Conv.ErrorBtwnSol2D(Results(1).Field, Results(2).Field, Results(3).Field, Results(4).Field, Results(5).Field,dx,dt,'Time')
table2latex(TField_T,'StraightSL00100100kd1_Field_Time')

[TLine_T]  = Conv.ErrorBtwnSol(Results(1).Line, Results(2).Line, Results(3).Line, Results(4).Line, Results(5).Line,drho,dt,'Time')
table2latex(TLine_T,'StraightSL00100100kd1_Line_Time')