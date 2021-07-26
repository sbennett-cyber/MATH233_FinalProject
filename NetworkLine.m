classdef NetworkLine
    % NETWORKLINE Is a class containing functions to solve a coupled
    % Network-Field system with up to three channels in the network. The
    % computational domain is assumed to be square. Verticies of the
    % network connecting to only one edge must lie on the boundary of the
    % computational domain. All other verticies must lie interior to the
    % computational domain. 
    %
    % The class has the following methods:
    %   NetworkLine -> Initializes the class and runs setup on the network.
    %
    %   SetIC -> Uses the properties of the class to set and return the initial conditions
    %
    %   PlotEBM     -> Plots the points used in the Embedded Boundary Method for easy debugging. 
    %
    %   CheckBoundary -> Identifies if a point lies on the boundary of the computational domain using properties of the class.
    %
    %   NormalLine -> Finds the intersection between a computational gridline in X or Y and a line. 
    %
    %   LagrIP_CGamma -> Uses a lagrange interpolating polynomial to find the coefficents for C_Gamma in the Embedded Boundary Method
    %
    %   LagrIP_CBar -> Uses a lagrange interpolating polynomial to find the coefficents for C_Bar in the Embedded Boundary Method
    %
    %
    %
    %
    %
    %
    %
    %
    %
    
    
    properties 
        Network;NetMtrix;
        Matrix; NumDomains; NumChannels
        minX;minY;maxX;maxY;Xm;Ym;BaseDrho;
        dx;dy;dt;tol;
        alpha; beta; CenterX; CenterY; IC; TrueDM5;
        Coeff_KPP1;Coeff_KPP2;KPP1;KPP2
    end
    
    methods 
        %% Initialize the Class
        function obj = NetworkLine(poly_, net_, BoundaryX, BoundaryY, PntsField, PntsLine, alpha_, beta_)
            
            obj.Network = net_; %Saves the initial structure of the network
            
            obj.NumDomains = length(poly_); %Saves the number of polygons
            obj.NumChannels = length(unique(net_.LineID)); %Saves the number of channels
            
            %Sets the size of the domain in X
            obj.minX = BoundaryX(1); 
            obj.maxX = BoundaryX(2);
            
            %Sets the size of the domain in Y
            obj.minY = BoundaryY(1); 
            obj.maxY = BoundaryY(2);
            
            %Create meshgrid for 2D problem
            X=linspace(BoundaryX(1),BoundaryX(2),PntsField(1)); 
            Y=linspace(BoundaryY(1),BoundaryY(2),PntsField(2));
            [XM,YM]=meshgrid(X,Y);
            obj.Xm = XM;
            obj.Ym = YM;
            
            %Set shape of gaussian IC
            obj.alpha = alpha_; 
            obj.beta = beta_;
            
            %Set center of gaussian IC
            obj.CenterX = net_.CenterGaussianX; 
            obj.CenterY = net_.CenterGaussianY;

            %Set function for IC
            obj.IC   = @(x,y) exp((-obj.alpha*(x-obj.CenterX).^2-obj.beta*(y-obj.CenterY).^2)); 
            %obj.IC   = @(x,y) 0*x+1; 
            obj.TrueDM5=@(x,y) exp(-obj.alpha*(x-obj.CenterX).^2-obj.beta*(y-obj.CenterY).^2).*((-2*obj.alpha)^2*(x-obj.CenterX).^2-(2*obj.alpha)+(-2*obj.beta)^2*(y-obj.CenterY).^2-(2*obj.beta)); %Set function for Laplacian
 
            %Find the spatial step size based on the size of the domain and the number of points
            obj.dx =abs(BoundaryX(2)-BoundaryX(1))/(PntsField(1)-1); 
            obj.dy =abs(BoundaryY(2)-BoundaryY(1))/(PntsField(2)-1);
            
            obj.BaseDrho = min(obj.dx,obj.dy)/2.25; %Set a basis for the spatial step on the network
            obj.tol = 10^-7; %Set tolerance
            
            obj.Coeff_KPP1 = poly_(1).r;
            obj.Coeff_KPP2 = poly_(2).r;
            
            obj.KPP1 = @(t,u) obj.Coeff_KPP1*u.*(1-u);
            obj.KPP2 = @(t,u) obj.Coeff_KPP2*u.*(1-u);
            
            % Setup the network in order to generate the matrix of coefficients later.
            for i = 1:obj.NumChannels %For each channel, find the length and discretize
                obj.NetMtrix(i).LineID = i; %Store the ID
                Index = find(obj.Network.LineID==i); %Find the index where the line is i
                vertexX = obj.Network.vertexX(Index(1),:); %Get the X vertex of the line
                vertexY = obj.Network.vertexY(Index(1),:); %Get the Y vertex of the line
                d = sqrt((vertexX(1)-vertexX(2))^2+(vertexY(1)-vertexY(2))^2); %Find the length of the line
                obj.NetMtrix(i).drho = d/(PntsLine(i)-1); %Correct drho based on the rounded number of points
                obj.NetMtrix(i).Size = PntsLine(i);
                obj.NetMtrix(i).Disc = linspace(0,d,PntsLine(i));
                m_and_b = polyfit(vertexX,vertexY,1);
                obj.NetMtrix(i).DiscIn2D_X(:,1) = linspace(vertexX(1),vertexX(2),PntsLine(i));
                obj.NetMtrix(i).DiscIn2D_Y(:,1) = m_and_b(1)*obj.NetMtrix(i).DiscIn2D_X(:)+m_and_b(2);  
                obj.NetMtrix(i).r = obj.Network.D(Index(1))/(2*obj.NetMtrix(i).drho^2); %Create the coefficient D/2drho. NOTE: This does not include time!
                obj.NetMtrix(i).D = obj.Network.D(Index(1)); %Create the coefficient D/2drho. NOTE: This does not include time!
                obj.NetMtrix(i).Sum_kd = sum(obj.Network.kd(Index)); %Find the domain that bounds this line
                obj.NetMtrix(i).Pointer_ka(:,:) = [obj.Network.Boundary(Index,1),obj.Network.ka(Index)']; %Saves the [domain, ka value]
                if ~isempty(obj.Network.SharedVertexX) %If there is a shared vertex
                    obj.NetMtrix(i).LocVertex = find(abs(obj.NetMtrix(i).DiscIn2D_X(:)- obj.Network.SharedVertexX)<obj.tol & abs(obj.NetMtrix(i).DiscIn2D_Y(:)- obj.Network.SharedVertexY)<obj.tol); %Identify which location has the shared vertex                 
                end
            end
        end
         
        %% Setup the initial conditions
        function [out1, out2] = GetXm(obj)
            out1 = obj.Xm;
            out2 = obj.Ym;
        end
        %% Setup the initial conditions
        function [out1, out2] = SetIC(obj,poly, netX, netY)
            %SETIC sets the initial conditions on the network and in each
            %subdomain. If there is more than one channel, it also returns
            %the IC at the junction. Returns two structures. 
            
            LineIC.Line(:,1) = obj.IC(netX,netY); %Find IC for all points on line
            if obj.NumChannels>1
                LineIC.Junction(:,:) = obj.IC(obj.Network.SharedVertexX,obj.Network.SharedVertexX);
            end

            for i = 1:obj.NumDomains
                FieldIC(i).Field(:,1) = obj.IC(poly(i).InXNB,poly(i).InYNB);
            end
            
            out1 = LineIC;
            out2 = FieldIC;
        end
        
        %% Plot solutions
        
        function [] = PlotSolutions(obj, poly, SolField)
            
            SizeX1=size(obj.Xm(2:end-1,2:end-1));
            SizeY1=size(obj.Ym(2:end-1,2:end-1));
            
            SizeX=size(obj.Xm);
            SizeY=size(obj.Ym);
            
         %   xm=reshape(obj.Xm',SizeX(1)*SizeX(2),1); %Use prime to ensure I reshape by row rather than column
         %   ym=reshape(obj.Ym',SizeY(1)*SizeY(2),1);
            u = sparse(SizeX1(1)*SizeX1(2),1);
            U = sparse(SizeX(1),SizeX(2));
            
            for i = 1:obj.NumDomains
                u(poly(i).LocateNB) = SolField(i).Field;
               % [in,bd] = inpolygon(poly(i).InXNB,poly(i).InYNB,poly(i).X,poly(i).Y);
               % PntsOnNet = find(bd==1);
               % u(poly(i).LocateNB(PntsOnNet)) = 0;
            end
            
            U(2:end-1,2:end-1) = reshape(u,SizeX1(1), SizeX1(2))';
            subplot(1,2,1)
                surf(obj.Xm, obj.Ym, U)
                view(2)
                axis equal
%             subplot(1,3,2)
%                 pcolor(obj.Xm, obj.Ym, U)
%                 %axis equal
            subplot(1,2,2)
                contour(obj.Xm,obj.Ym, U,[0.01,0.01],'LineWidth',2)
                axis equal
            
            
        end
        
        
        %% Plot the results of the EBM
        function [] = PlotEBM(obj,poly, CenterID, Center, Neigbors, C_Gamma, C_i, Interp)
            % PLOTEBM Plots the points identified by the Embedded Boundary
            % Method. This function is automiatically triggered anytime the
            % EBM is used in a function in the NetworkLine Class.
            % 
            % Subdomians 1, 2, and 3 are plotted in red, blue, and green
            % with the corresponding discretization points in each domain
            % plotted in the same color. Note that points shared by more
            % than one domain may not plot in the correct color. 
            %
            % The center of the five point stencil is denoted by a back
            % circle. 
            %
            % The point in the stencil falling outside the domain
            % of interest is a black astrisk. 
            %
            % A white astrisk denotes the intersection of the closest
            % boundary and the line segment normal to the boundary passing
            % through the missing stencil point. 
            %
            % Yellow astrisks indicate where the normal passes through
            % either either gridlines in X or Y within the domain of
            % interest.
            %
            % Magenta points identify gridpoints used to get the values at
            % the yellow points. 
            figure
               for i = 1:obj.NumDomains
                    if i == 1
                        plot(poly(1).pgon,'FaceColor','green')
                        hold on
                        plot(poly(1).InX,poly(1).InY, 'g*')
                    elseif i == 2
                        plot(poly(2).pgon,'FaceColor','blue')
                        hold on
                        plot(poly(2).InX,poly(2).InY, 'b*')
                    elseif i == 3
                        plot(poly(3).pgon,'FaceColor','red')
                        hold on
                        plot(poly(3).InX,poly(3).InY, 'r*')
                    else
                        "More than three subdomains"
                        plot(obj.poly(4).pgon,'FaceColor','cyan')
                        hold on
                        plot(poly(4).InX,poly(4).InY, 'c*')
                    end
               end
               
               plot(Center(1), Center(2),'ko')
               
               if ~isempty(Neigbors)
                    plot(Neigbors(:,1), Neigbors(:,2),'k*')
               end
               if ~isempty(C_Gamma)
                    plot(C_Gamma(:,1), C_Gamma(:,2),'w*')
               end
               if ~isempty(C_i)
                   plot(C_i(:,1), C_i(:,2),'c*')
               end
               if ~isempty(Interp)
                   plot(Interp(:,1), Interp(:,2),'m*')
               end
               
               axis equal
                xlabel("X")
                ylabel("Y")
                drawnow
                hold off
        end
        
        %% Check if a point lies on the boundary
        function [XOnBound, YOnBound] = CheckBoundary(obj,X,Y)
            %CHECKBOUNDARY returns 1 if the input point lies on the
            %boundary of the computational domain as defined by the
            %properties of the class.
            XOnBound = find(abs(X - obj.minX)<obj.tol | abs(X - obj.maxX)<obj.tol);    %Find the location of points on the boundary
            YOnBound = find(abs(Y - obj.minY)<obj.tol | abs(Y - obj.maxY)<obj.tol);
        end
        
        %% Given a X or Y gridline, find where the normal intersects the gridline
        function [out] = NormalLine(obj,m,b,u,Case)
            %NORMALLINE takes in a slope m, y-intercept b, and either x or
            %y part of a point (x,y) and returns either y or x to complete the
            %coordinate. 
            switch Case
                case "x"
                    out=((u-b)/m);
                case "y"
                    out=m*u+b;
            end
        end
        
        %% Coefficents for C_Gamma
        function [out] = LagrIP_CGamma(obj, YSize, xr, x, ka, kd, D)
            % LAGRIP_CGAMMA uses Lagrange Interpolating Polynomials to calculate
            %the coefficents when solving for C_Gamma in the embedded
            %boundary method. 

            switch YSize 
                case 1
                    A_omega = 0;
                    B_omega = 0;
                    C_omega = 0;
                    D_omega = (1/xr);
                    Lower = D*D_omega+ka;
                case 2
                    A_omega = (xr+x)/(x*xr);
                    B_omega = (xr)/(x^2+x*xr);
                    C_omega = 0;
                    D_omega = (1/xr) + (1/(xr+x));
                    Lower = D*D_omega+ka;
                case 3
                    A_omega = (xr^2+xr*x+2*xr*x+2*x^2)/(2*x^2*xr); %Same as Shilpa's
                    B_omega = (xr^2+2*x*xr)/(x^3+x^2*xr);          %Same as Shilpa's
                    C_omega = (xr^2+  x*xr)/(4*x^3+2*x^2*xr);      %Same as Shilpa's
                    D_omega = (1/xr) + (1/(xr+x)) + (1/(xr+2*x));  %Different!
                    Lower = D*D_omega+ka;               
            end
            try
            out=[D*A_omega/Lower,-D*B_omega/Lower,D*C_omega/Lower,kd/Lower]; %Coeff on C_i, C_ii, C_iii, rho
            catch
                "Input too long"
                pause
            end
        end

        %% Coefficents for C_Bar
        function [out] = LagrIP_CBar(obj, YSize, xbar, xr, x, Coeff_CGamma)
            %LAGRIP_CBAR uses Lagrange Interpolating Polynomials to calculate
            %the coefficents when solving for C_Gamma in the embedded boundary method. 
            switch YSize
                case 1
                    A = -xbar/xr;
                    B = 0;
                    C = 0;
                    D = x/xr;
                case 2
                    A = (-2*xbar)/(xr);
                    B = (xbar)/(x+xr);
                    C = 0;
                    D = (2*x^2)/((xr^2+xr*x));
                case 3
                    A = (-3*xbar)/(xr);
                    B = (3*xbar)/(x+xr);
                    C = (-xbar)/(2*x+xr);
                    D = (6*x^3)/((xr^2+xr*x)*(xr+2*x));
            end
            out=[A+D*Coeff_CGamma(1),B+D*Coeff_CGamma(2),C+D*Coeff_CGamma(3),D*Coeff_CGamma(end)]; %Coeff on C_i, C_ii, C_ii, C_Gamma
        end
        
        %% Regular Lagrange Interpolating Polynomial
        function [out] = LagrangeIP(obj, Gridpnts, Crosspnt)
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here

            Case=length(Gridpnts);
            

                switch Case
                    case 1
                        out=1;
                    case 2
                        Coeff1=(Crosspnt-Gridpnts(2))/(Gridpnts(1)-Gridpnts(2));
                        Coeff2=(Crosspnt-Gridpnts(1))/(Gridpnts(2)-Gridpnts(1));
                        out=[Coeff1,Coeff2];
                    case 3
                        Coeff1=((Crosspnt-Gridpnts(2))*(Crosspnt-Gridpnts(3)))/((Gridpnts(1)-Gridpnts(2))*(Gridpnts(1)-Gridpnts(3)));
                        Coeff2=((Crosspnt-Gridpnts(1))*(Crosspnt-Gridpnts(3)))/((Gridpnts(2)-Gridpnts(1))*(Gridpnts(2)-Gridpnts(3)));
                        Coeff3=((Crosspnt-Gridpnts(1))*(Crosspnt-Gridpnts(2)))/((Gridpnts(3)-Gridpnts(1))*(Gridpnts(3)-Gridpnts(2)));
                        out=[Coeff1,Coeff2,Coeff3];
                    case 4
                        Coeff1=((Crosspnt-Gridpnts(2))*(Crosspnt-Gridpnts(3))*(Crosspnt-Gridpnts(4)))/((Gridpnts(1)-Gridpnts(2))*(Gridpnts(1)-Gridpnts(3))*(Gridpnts(1)-Gridpnts(4)));
                        Coeff2=((Crosspnt-Gridpnts(1))*(Crosspnt-Gridpnts(3))*(Crosspnt-Gridpnts(4)))/((Gridpnts(2)-Gridpnts(1))*(Gridpnts(2)-Gridpnts(3))*(Gridpnts(2)-Gridpnts(4)));
                        Coeff3=((Crosspnt-Gridpnts(1))*(Crosspnt-Gridpnts(2))*(Crosspnt-Gridpnts(4)))/((Gridpnts(3)-Gridpnts(1))*(Gridpnts(3)-Gridpnts(2))*(Gridpnts(3)-Gridpnts(4)));
                        Coeff4=((Crosspnt-Gridpnts(1))*(Crosspnt-Gridpnts(2))*(Crosspnt-Gridpnts(3)))/((Gridpnts(4)-Gridpnts(1))*(Gridpnts(4)-Gridpnts(2))*(Gridpnts(4)-Gridpnts(3)));
                        out=[Coeff1,Coeff2,Coeff3,Coeff4]; 
                end
                
        end
        
         %% Find the correct set of X and Y coordinates to use
        function [out1, out2, out3] = GenerateStencil(obj, Center, Cross, C_Gamma, Type, m)
            % Note here that m is the slope of the normal line
            % out3 sets the flux in the normal direction to be negative if y is decreasing 
            
            CheckX = Center(1)-Cross(1);
            CheckY = Center(2)-Cross(2);
            
            if (abs(Center(1)-C_Gamma(1))<obj.tol && abs(Center(2)-C_Gamma(2)) < obj.tol) %Check if Center equals Cross
                CenEquCGamm = 1;
            else 
                CenEquCGamm = 0;
            end
            
            if (abs(C_Gamma(1)-Cross(1))<obj.tol && abs(C_Gamma(2)-Cross(2)) < obj.tol) %Check if C_Gamma equals Cross
                CGammEquCross = 1;
            else 
                CGammEquCross = 0;
            end
            
            
            
            X = obj.Xm(1,:)';
            Y = obj.Ym(:,1);
            
            if (abs(Center(2)-C_Gamma(2))<obj.tol && abs(Center(1)-C_Gamma(1))<obj.tol && abs(Center(2)-Cross(2))>obj.tol)  %Check if Center and C_Gamma are the same
                out3 = "Equal"; %If they are the same, return equal
            else
                out3 = "NotEqual"; %If they are not, return not equal
            end
            
            LocXTemp = find((abs(X-Cross(1))-min(abs(X-Cross(1))))<obj.tol); %Find the location of the x coordinate of Cross in X discretization
            LocYTemp = find((abs(Y-Cross(2))-min(abs(Y-Cross(2))))<obj.tol); %Find the location of the y coordinate of Cross in Y discretization
                
            if (X(LocXTemp(1))-Cross(1))>obj.tol %If X at LocX is greater than Cross, shift to the left
                LocX = X(LocXTemp(1)-1);
            else 
                LocX = X(LocXTemp(1));
            end
            if (Y(LocYTemp(1))-Cross(2))>obj.tol %If Y at LocY is greater than Cross, shift to the left
                LocY = Y(LocYTemp(1)-1);
            else 
                LocY = Y(LocYTemp(1));
            end
                    
        
                switch Type
                    case "Horizontal"
                        x_iStencil = [];
                        if CGammEquCross == 1 && CenEquCGamm == 1
                            y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY-obj.dy, LocY-2*obj.dy,LocY-3*obj.dy];

                        elseif CGammEquCross == 0 && CenEquCGamm == 1
%                             if CheckY > 0
%                                 y_iStencil = [LocY+2*obj.dy,LocY+obj.dy,LocY];
%                             else
%                                 y_iStencil = [LocY,LocY-obj.dy,LocY-2*obj.dy];
%                             end
%                             if CheckY > 0
%                                 y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy];
%                             else
%                                 y_iStencil = [LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
%                             end
                            if CheckY > 0
                                y_iStencil = [LocY+4*obj.dy,LocY+3*obj.dy,LocY+2*obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy];
                            else
                                y_iStencil = [LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-2*obj.dy,LocY-3*obj.dy,LocY-4*obj.dy];
                            end

                            
                        else
                            if CheckY > 0
                                y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy];
                            else
                                y_iStencil = [LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
                            end

                        end
                    case "Angled"
                        
                        if CGammEquCross == 1 && CenEquCGamm == 1
                            if abs(X(LocXTemp(1))-Cross(1))<obj.tol
                                x_iStencil = [LocX+3*obj.dx,LocX+2*obj.dx,LocX+obj.dx,LocX-obj.dx,LocX-2*obj.dx,LocX-3*obj.dx];
                               % y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
                            
                            else
                                x_iStencil = [LocX+3*obj.dx,LocX+2*obj.dx,LocX+obj.dx,LocX,LocX-obj.dx,LocX-2*obj.dx]; %Changed July 14, 2021
                               % y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy];
                            end
                            if abs(Y(LocYTemp(1))-Cross(2))<obj.tol
                               % x_iStencil = [LocX+3*obj.dx,LocX+2*obj.dx,LocX+obj.dx,LocX-obj.dx,LocX-2*obj.dx,LocX-3*obj.dx];
                                y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
                            
                            else
                               % x_iStencil = [LocX+3*obj.dx,LocX+2*obj.dx,LocX+obj.dx,LocX,LocX-obj.dx,LocX-2*obj.dx]; %Changed July 14, 2021
                                y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy];
                            end
                            
                            
                        else
                            if CheckY > obj.tol
                                y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy];
                            elseif abs(CheckY) < obj.tol
                                y_iStencil = [LocY+3*obj.dy,LocY+2*obj.dy,LocY+obj.dy,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
                            else
                                y_iStencil = [LocY+2*obj.dy,LocY+obj.dy,LocY,LocY-obj.dy,LocY-2*obj.dy,LocY-3*obj.dy];
                            end
                            
                            if CheckX > obj.tol
                                x_iStencil = [LocX+3*obj.dx,LocX+2*obj.dx,LocX+obj.dx,LocX,LocX-obj.dx,LocX-2*obj.dx];
                            elseif abs(CheckX) < obj.tol
                                x_iStencil = [LocX+3*obj.dy,LocX+2*obj.dx,LocX+obj.dx,LocX-obj.dx,LocX-2*obj.dx,LocX-3*obj.dx];
                            else
                                x_iStencil = [LocX+2*obj.dx,LocX+obj.dx,LocX,LocX-obj.dx,LocX-2*obj.dx,LocX-3*obj.dx];
                            end

                        end
                end
                        
            out1 = x_iStencil;
            out2 = y_iStencil;
            
            
        end
        %% Find the correct set of X and Y coordinates to use
        function [out1, out2, out3] = GenerateStencil_Old(obj, Center, Cross, C_Gamma, Type, m)
            % Note here that m is the slope of the normal line
            % out3 sets the flux in the normal direction to be negative if y is decreasing 
            
            CheckX = Center(1)-Cross(1);
            CheckY = Center(2)-Cross(2);
            
            X = obj.Xm(1,:)';
            Y = obj.Ym(:,1);
            
            if (abs(Center(2)-Cross(2))<obj.tol && abs(Center(1)-Cross(1))<obj.tol) && (abs(Center(2)-C_Gamma(2))<obj.tol && abs(Center(1)-C_Gamma(1))<obj.tol) %Check if Center and Cross are the same
                out3 = "NotEqual"; %If they are not, return not equal
            elseif (abs(Center(2)-C_Gamma(2))<obj.tol && abs(Center(1)-C_Gamma(1))<obj.tol) && (abs(Center(2)-Cross(2))>obj.tol && abs(Center(1)-Cross(1))>obj.tol)  %Check if Center and C_Gamma are the same
                out3 = "Equal"; %If they are the same, return equal
            elseif (abs(Center(2)-C_Gamma(2))<obj.tol && abs(Center(1)-C_Gamma(1))<obj.tol) && (abs(Center(2)-Cross(2))>obj.tol && abs(Center(1)-Cross(1))<obj.tol) && Type == "Horizontal"  %Check if Center and C_Gamma are the same
                out3 = "Equal"; %If they are the same, return equal
            else
                out3 = "NotEqual"; %If they are not, return not equal
            end
            
            if min(abs(X-Cross(1)))<obj.tol && min(abs(Y-Cross(2)))<obj.tol %Check of Cross lies on a gridpoint
                FixedPnt = Cross;
        
                switch Type
                    case "Horizontal"
                        x_iStencil = [];
                        if CheckY>0
                            y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                        elseif abs(CheckY)<obj.tol
                            y_iStencil = [FixedPnt(2)+3*obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy, FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];
                            
                        else
                            y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                        end
                    case "Angled"
                        if abs(CheckX)<obj.tol && abs(CheckY)<obj.tol 
                            x_iStencil = [FixedPnt(1)+3*obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+obj.dx, FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx,FixedPnt(1)-3*obj.dx]; %New June 14
                            y_iStencil = [FixedPnt(2)+3*obj.dx, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy, FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy,FixedPnt(2)-3*obj.dy];
                        elseif abs(CheckX)<obj.tol 
                            if m > 0 && CheckY<0    
                               x_iStencil = [FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx, FixedPnt(1)-3*obj.dx];
                               y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy]; 

                            elseif m > 0 && CheckY>0
                               x_iStencil = [FixedPnt(1)+obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+3*obj.dx];
                               y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                            elseif m < 0 && CheckY<0
                               x_iStencil = [FixedPnt(1)+obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+3*obj.dx];
                               y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                            elseif m < 0 && CheckY>0
                               x_iStencil = [FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx, FixedPnt(1)-3*obj.dx];
                               y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];
                            else
                                "Error in Generate Stencil. Case in CheckX==0 does not exist."
                            end

                        elseif abs(CheckY)<obj.tol
                            if m > 0 && CheckX<0    
                               x_iStencil = [FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx, FixedPnt(1)-3*obj.dx];
                               y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, Cross(2)-3*obj.dy]; 

                            elseif m > 0 && CheckX>0
                               x_iStencil = [FixedPnt(1)+obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+3*obj.dx];
                               y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                            elseif m < 0 && CheckX>0
                               x_iStencil = [FixedPnt(1)+obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+3*obj.dx];
                               y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                            elseif m < 0 && CheckX<0
                               x_iStencil = [FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx, FixedPnt(1)-3*obj.dx];
                               y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];
                            else
                                "Error in Generate Stencil. Case in CheckX==0 does not exist."
                            end
                        else
                            "Error in Generate Stencil"
                        end
                end
                        
            elseif min(abs(X-Cross(1)))<obj.tol %Check if it lies on a gridline in X
        %%%        "On X Gridline"
                 LocMinX = find((abs(X-Cross(1)) - min(abs(X-Cross(1))))<obj.tol & X>obj.tol); %Find the location of the closest gridline in X
                 FixedPnt(1) = X(LocMinX(1));
                 
                 LocMinY = find((abs(Y-Cross(2)) - min(abs(Y-Cross(2))))<obj.tol & Y>obj.tol); %Find the location of the closest gridline in Y      
                 if Y(LocMinY(1))>Cross(2)
                    FixedPnt(2) = Y(LocMinY(1))-obj.dy;
                else
                    FixedPnt(2) = Y(LocMinY(1));
                end

                 switch Type
                    case "Horizontal"
                        x_iStencil = [0,0,0];
                        if CheckY>0
                            y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                        elseif abs(CheckY)<obj.tol
                            y_iStencil = [FixedPnt(2)+3*obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy, FixedPnt(2), FixedPnt(2)-1*obj.dy, FixedPnt(2)-2*obj.dy];
                            
                        else
                            y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                        end
                    case "Angled"                     

                        x_iStencil = [FixedPnt(1)+3*obj.dx, FixedPnt(1)+2*obj.dx, FixedPnt(1)+obj.dx, FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx,FixedPnt(1)-3*obj.dx];
                        y_iStencil = [FixedPnt(2)+3*obj.dx, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy,   FixedPnt(2),        FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy];
                 end
                 
             elseif min(abs(Y-Cross(2)))<obj.tol %Check if it lies on a gridline in Y
     %%%           "On Y Gridline"
                LocMinY = find((abs(Y-Cross(2)) - min(abs(Y-Cross(2))))<obj.tol & Y>obj.tol); %Find the location of the closest gridline in Y
                FixedPnt(2) = Y(LocMinY(1));
                 switch Type
                    case "Horizontal"
                        x_iStencil = [0,0,0];
                        if CheckY>0
                            y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                        elseif abs(CheckY)<obj.tol
                            y_iStencil = [FixedPnt(2)+3*obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy, FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];
                            
                        else
                            y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                        end
                    case "Angled"
                        LocMinX = find((abs(X-Cross(1)) - min(abs(X-Cross(1))))<obj.tol & X>obj.tol); %Find the location of the closest gridline in X
                        if X(LocMinX(1))>Cross(1)
                            FixedPnt(1) = X(LocMinX(1))-obj.dx;
                        else 
                            FixedPnt(1) = X(LocMinX(1));
                        end

                        y_iStencil = [FixedPnt(2)+3*obj.dx, FixedPnt(2)+2*obj.dx, FixedPnt(2)+obj.dx, FixedPnt(2)-obj.dx, FixedPnt(2)-2*obj.dx, FixedPnt(2)-3*obj.dx];
                        x_iStencil = [FixedPnt(1)+3*obj.dy, FixedPnt(1)+2*obj.dy, FixedPnt(1)+obj.dy,   FixedPnt(1),        FixedPnt(1)-obj.dy, FixedPnt(1)-2*obj.dy];
                 end
            else
     %%%           "Not On a Gridline"
                LocMinY = find((abs(Y-Cross(2)) - min(abs(Y-Cross(2))))<obj.tol & Y>obj.tol); %Find the location of the closest gridline in Y

                if Y(LocMinY(1))>Cross(2)
                    FixedPnt(2) = Y(LocMinY(1))-obj.dy;
                else
                    FixedPnt(2) = Y(LocMinY(1));
                end
                
                 switch Type
                    case "Horizontal"
                        x_iStencil = [0,0,0];
                        if CheckY>0
                            y_iStencil = [FixedPnt(2)+obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+3*obj.dy];

                        elseif abs(CheckY)<obj.tol
                            y_iStencil = [FixedPnt(2)+3*obj.dy, FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy, FixedPnt(2), FixedPnt(2)-1*obj.dy, FixedPnt(2)-2*obj.dy];
                            
                        else
                            y_iStencil = [FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy, FixedPnt(2)-3*obj.dy];

                        end
                    case "Angled"
                       LocMinX = find((abs(X-Cross(1)) - min(abs(X-Cross(1))))<obj.tol & X>obj.tol); %Find the location of the closest gridline in X
                       

                       if X(LocMinX(1))>Cross(1)
                            FixedPnt(1) = X(LocMinX(1))-obj.dx;
                       else 
                            FixedPnt(1) = X(LocMinX(1));
                       end
                       x_iStencil = [FixedPnt(1)+3*obj.dx,FixedPnt(1)+2*obj.dx, FixedPnt(1)+obj.dx,   FixedPnt(1),        FixedPnt(1)-obj.dx, FixedPnt(1)-2*obj.dx];
                       y_iStencil = [FixedPnt(2)+3*obj.dy,FixedPnt(2)+2*obj.dy, FixedPnt(2)+obj.dy,   FixedPnt(2),        FixedPnt(2)-obj.dy, FixedPnt(2)-2*obj.dy];
                 end
            end
            
            
            out1 = x_iStencil;
            out2 = y_iStencil;
            
            
        end
        
        %% Find if the line falls between Center and Crosspoint
        function [out] = FindChannel(obj, v1, v2, Center, Cross)
              if abs(v1(2)-v2(2))<obj.tol                       %Check if line segment is a horizontal line
                  out = ((v1(2)> min(Center(2), Cross(2)) && v1(2)< max(Center(2), Cross(2))) || abs(min(Center(2), Cross(2))-v1(2))<obj.tol ||abs(max(Center(2), Cross(2))-v1(2))<obj.tol); %Check that the y of the line falls between Center and Cross;;;
              elseif abs(v1(1)-v2(1))<obj.tol                   %Check if vertical line segment
                  "ERROR! Boundary cannot be vertical"
              else
                  LineSpecs = polyfit([v1(1),v2(1)],[v1(2),v2(2)],1); %Returns [m,b] for line
                  if abs(Center(1)-Cross(1))<obj.tol
                      y=Center(1)*LineSpecs(1)+LineSpecs(2);
                      out = ((y> min(Center(2), Cross(2)) && y< max(Center(2), Cross(2))) || abs(min(Center(2), Cross(2))-y)<obj.tol ||abs(max(Center(2), Cross(2))-y)<obj.tol); %Check that the y of the line falls between Center and Cross;;
                      
                  elseif abs(Center(2)-Cross(2))<obj.tol
                      x=(Center(2)-LineSpecs(2))/LineSpecs(1);
                      out = ((x> min(Center(1), Cross(1)) && x< max(Center(1), Cross(1))) || abs(min(Center(1), Cross(1))-x)<obj.tol ||abs(max(Center(1), Cross(1))-x)<obj.tol); %Check that the y of the line falls between Center and Cross;
                  end
              end
        end
        
        %% Find the normal to a channel of the network and passing through a crosspoint
        function [d, Intercept, Type, m_and_b] = point_to_line(obj, pt, v1, v2, Case)
            %%% This function takes in a point (pt) and the endpoints of a line segment
            %%% (v1,v2), and returns the shortest distance from the point to line along
            %%% the normal as well as the intercept of the normal and line segment. If
            %%% the intercept of the normal occcurs outside of the line segment, then 
            %%% the function will return the distance to the closest end point. 

            %%% NOTE: This function assumes that v1(1)<v2(1)
                  if max(size(v1))==2
                      v1 = [v1(1),v1(2),0];
                      v2 = [v2(1),v2(2),0];
                  end
                  if max(size(pt))==2
                      pt = [pt(1),pt(2),0];
                  end
                  a = v1 - v2;
                  b = pt - v2;
                  d = norm(cross(a,b)) / norm(a); %Find the shortest distance from point to line on the normal

                  if abs(v1(2)-v2(2))<obj.tol                       %Check if line segment is a horizontal line
                      if abs(v1(2)-pt(2))>obj.tol && v1(2)>pt(2)                    %If yes, the intercept can be found easily
                        Intercept = [pt(1),pt(2)+d];
                      elseif abs(v1(2)-pt(2))>obj.tol && v1(2)<pt(2)
                          Intercept = [pt(1),pt(2)-d];
                      elseif abs(v1(2)-pt(2))<obj.tol
                          Intercept = [pt(1),pt(2)];
                      end

                  elseif abs(v1(1)-v2(1))<obj.tol                   %Check if vertical line segment
                      "ERROR! Boundary cannot be vertical"
                  else
                      LineSpecs = polyfit([v1(1),v2(1)],[v1(2),v2(2)],1); %Returns [m,b] for line
                      a1=LineSpecs(1);
                      b1=LineSpecs(2);
                      a2 = -1/a1;
                      b2 = (1/a1)*pt(1)+pt(2);
                      Intercept = [(-b2+b1)/(a2-a1),a1*((-b2+b1)/(a2-a1))+b1];
                  end

                  if v1(1)<v2(1)                            %Identify which point has the largest and smallest coordinate in x
                      MinPnt = v1;
                      MaxPnt = v2;
                  else
                      MinPnt = v1;
                      MaxPnt = v2;
                  end

                  if Intercept(1)<MinPnt(1)           %Check if the intercept falls to the left of v1 in the X direction
                      Intercept = MinPnt;
                      d = sqrt((pt(1)-MinPnt(1))^2 + (pt(2)-MinPnt(2))^2);
                  elseif Intercept(1)>MaxPnt(1)       %Check if the intercept falls to the right of v2 in the X direction
                      Intercept = MaxPnt;
                      d = sqrt((pt(1)-MaxPnt(1))^2 + (pt(2)-MaxPnt(2))^2);
                  end   

                  if abs(v1(2)-v2(2))<obj.tol %Check if the boundary is horizontal. This always means the normal is vertical
                      Type = "Horizontal";
                      m_and_b = [NaN,NaN];
                  
                  elseif abs(v1(2)-v2(2))>obj.tol && abs(Intercept(1)-pt(1))<obj.tol && abs(d)>obj.tol%Check if the normal is vertical even if the line is angled. 
                      Type = "Horizontal";
                      m_and_b = [NaN,NaN];
                      
                  elseif abs(v1(2)-v2(2))>obj.tol && abs(d)<obj.tol %Check if line segment is angled and pt falls on line segment
                      Type = "Angled";
                      m_and_b = [a2,b2];
                  else                                      %Check if line segment is angled and pt does not fall on line segment
                      Type = "Angled";
                      m_and_b = polyfit([Intercept(1),pt(1)],[Intercept(2),pt(2)],1); %Returns [m,b] for line
                  end

                  switch Case
                      case "Plot"
                          plot(pt(1),pt(2), 'r*')
                          hold on
                          plot([v1(1),v2(1)],[v1(2),v2(2)], 'Linewidth',2)
                          hold on
                          plot(Intercept(1),Intercept(2),'b*')
                          hold on
                          plot([pt(1),Intercept(1)],[pt(2),Intercept(2)], 'Linewidth',2)
                          xlabel("X")
                          ylabel("Y")
                          legend ("Point of Interest", "Initial Line Segment", "Intercept between Initial Line Segment and Normal","Line Segment along Normal")
                      case "No Plot"
                  end


        end
            
                %% Find the location of C_i, C_ii, and C_iii and return the coordinates & coefficients
        function [out1, out2] = Find_CBar(obj, xr_bar, C_Gamma, X_Cross, XSize, poly, CenterID, ka, kd, D, Case)
            %%% NOTE: There is a flaw in this function. If the value C_i,
            %%% C_ii, or C_iii lies on a gridpoint, I do not test for this.
            %%% Instead, I interpolate and if I do not have enough points I
            %%% set to zero. This can be a problem in some cases and should
            %%% be corrected. May 29, 2021.
            
            
            if XSize(1)==0                                                         % if there are no points in X or only one point in X, do nothing
                CoeffRho_CGamma = 0;
                CoeffRho_CBar   = 0;
                MtrixCorrection_CBar = [0,0,0];
                
            else 
                if sqrt(abs(C_Gamma(1,1)-X_Cross(1,1)).^2+abs(C_Gamma(1,2)-X_Cross(1,2)).^2)>sqrt(abs(C_Gamma(1,1)-X_Cross(end,1)).^2+abs(C_Gamma(1,2)-X_Cross(end,2)).^2)
                    X_Cross = flip(X_Cross);
                end
                
                MtrixCorrection_CBar = [];
                xi_r = sqrt(min(abs(C_Gamma(1,1)-X_Cross(:,1))).^2+min(abs(C_Gamma(1,2)-X_Cross(:,2))).^2);       % Find the distance between x_i and x_r -> Fixed on June 4
                if XSize(1)==1
                    xi_ii = 0;
                else 
                    xi_ii   = sqrt(abs(X_Cross(2,1)-X_Cross(1,1)).^2+abs(X_Cross(2,2)-X_Cross(1,2)).^2);                                                                    % Find the distance between x_i and x_r
                end
                
                Coeff_CGamma    = obj.LagrIP_CGamma(XSize(1), xi_r, xi_ii, ka, kd, D);                    % Calculate Lagrange IP coefficents for C_Gamma
                CoeffRho_CGamma = Coeff_CGamma(end);                                              % Identify the coefficent on the rho term
                
%                 if abs(sum(Coeff_CGamma)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                       "Coeff_Ci Greater than 1"
%                       pause
%                 end
                
                if abs(xr_bar)>obj.tol
                    Coeff_CBar = obj.LagrIP_CBar(XSize(1), xr_bar, xi_r, xi_ii, Coeff_CGamma); % Calculate Lagrange IP coefficents for C_bar
                else 
                    Coeff_CBar = Coeff_CGamma;                                                  % If C_bar = C_Gamma, use the same Lagrange IP
                end 
                CoeffRho_CBar = Coeff_CBar(end); %Identify coefficient on rho term
%                 if abs(sum(Coeff_CBar)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                       "Coeff_Ci Greater than 1"
%                      % pause
%                 end
                
                switch Case
                    case "X"
                        Y_Below1=X_Cross(:,2)-mod(X_Cross(:,2),obj.dx);                                          % Find the gridpoint below or to the left of each cross point
                        Y_Vals=[Y_Below1-2*obj.dx,Y_Below1-obj.dx,Y_Below1,Y_Below1+obj.dx,Y_Below1+2*obj.dx,Y_Below1+3*obj.dx]; % Find six points in x surrounding each point in y 
                        X_Vals = Y_Vals*0 + X_Cross(:,1);
                        Y_ValsExtraL = Y_Below1-3*obj.dx;
                        Y_ValsExtraR = Y_Below1+4*obj.dx;
                        X_ValsExtraL = X_Cross(:,1);
                        X_ValsExtraR = X_Cross(:,1);
                        
                        
                    case "Y"
                        Y_Below1=X_Cross(:,1)-mod(X_Cross(:,1),obj.dx);                                          % Find the gridpoint below or to the left of each cross point
                        X_Vals=[Y_Below1-2*obj.dx,Y_Below1-obj.dx,Y_Below1,Y_Below1+obj.dx,Y_Below1+2*obj.dx,Y_Below1+3*obj.dx]; % Find six points in x surrounding each point in y
                        Y_Vals = X_Vals*0 + X_Cross(:,2);
                        X_ValsExtraL = Y_Below1-3*obj.dx;
                        X_ValsExtraR = Y_Below1+4*obj.dx;
                        Y_ValsExtraL = X_Cross(:,2);
                        Y_ValsExtraR = X_Cross(:,2);
                end

                for j=1:XSize(1)
                    
                    [Ind, ~] = inpolygon(X_Vals(j,:),Y_Vals(j,:),poly(CenterID).X,poly(CenterID).Y); %Identify which points fall within the boundary
                    
                    FinalValues = [X_Vals(j,:)',Y_Vals(j,:)',zeros(length(Y_Vals(j,:)),1)]; %Store values as a matrix of form [X,Y]
                    FinalValues=FinalValues(Ind'==1 & FinalValues(:,1)>obj.minX & FinalValues(:,1)<obj.maxX & FinalValues(:,2)>obj.minY & FinalValues(:,2)<obj.maxY,:);  %Remove all points that fall outside the domain of interest or on the boundary
                    
                    ExtraL = [X_ValsExtraL(j,:),Y_ValsExtraL(j,:)];
                    [IndL, ~] = inpolygon(X_ValsExtraL(j,:),Y_ValsExtraL(j,:),poly(CenterID).X,poly(CenterID).Y); %Identify which points fall within the boundary
                    ExtraL=ExtraL(IndL'==1 & ExtraL(:,1)>obj.minX & ExtraL(:,1)<obj.maxX & ExtraL(:,2)>obj.minY & ExtraL(:,2)<obj.maxY,:);
                    
                    ExtraR = [X_ValsExtraR(j,:),Y_ValsExtraR(j,:)]; 
                    [IndR, ~] = inpolygon(X_ValsExtraR(j,:),Y_ValsExtraR(j,:),poly(CenterID).X,poly(CenterID).Y); %Identify which points fall within the boundary
                    ExtraR=ExtraR(IndR'==1 & ExtraR(:,1)>obj.minX & ExtraR(:,1)<obj.maxX & ExtraR(:,2)>obj.minY & ExtraR(:,2)<obj.maxY,:);
                    
                    switch Case
                        case "X"
                            CheckDist = find((abs(FinalValues(1:end-1,2)-FinalValues(2:end,2))-obj.dx)>obj.tol);
                        case "Y"
                            CheckDist = find((abs(FinalValues(1:end-1,1)-FinalValues(2:end,1))-obj.dx)>obj.tol);
                    end
                    
                    
                    if ~isempty(CheckDist) & length(FinalValues)>1 %Check to see if a point falls into another domain
                       % Cut = ones(length(FinalValues),1);
                        for q = 1:length(CheckDist)
%                             switch Case
%                                 case "X"
%                                     First = abs(FinalValues(CheckDist(q),2)-X_Cross(j,2));  %Find the distance between the crosspoing and the endpoints of FinalValues
%                                     Last =  abs(FinalValues(CheckDist(q)+1,2)-X_Cross(j,2));
%                                 case "Y"
%                                     First = abs(FinalValues(CheckDist(q),1)-X_Cross(j,1));
%                                     Last =  abs(FinalValues(CheckDist(q)+1,1)-X_Cross(j,1));
%                             end

                           % LocMax = find(abs(Cut-max(Diff))<obj.tol);
                            %LocMin = find(abs(Cut-min(Diff))<obj.tol);
                            Location = find(Ind==0);
                            if sum(Ind(1:Location(1))) > sum(Ind(Location(1):end))
                                FinalValues = FinalValues(1:Location-1,:);
                                Side = 'L';
                            elseif sum(Ind(1:Location(1))) < sum(Ind(Location(1):end)) %#ok<*AND2>
                                FinalValues = FinalValues(Location+1:end,:);
                                Side = 'R';
                            else
                                if length(Y_Cross(:,1))==6
                                    FinalValues = FinalValues(1:3,:);
                                end
                               % "XSize Incorrect" %#ok<*NOPRT>
                                %pause
                            %pause
                            end
%                             if First>Last %If first is farther away from the point we want to interpolate
%                                 Cut(1:CheckDist(q)) = 0; %Remove the first point
%                                 Side = 'L';
%                             else
%                                 Cut(CheckDist(q)+1:end) = 0; %Remove the last point
%                                 Side = 'R';
%                             end   
                            
                        end
                        %FinalValues = FinalValues(Cut==1,:);
                        if length(FinalValues)<4
                            switch Side
                                case 'R'
                                    FinalValues = [FinalValues; [ExtraR,0]];
                                case 'L'
                                    FinalValues = [[ExtraL,0]; FinalValues];
                            end
                        end
                    end
                    
                    
                    if length(FinalValues(:,1)) == 6 %If there are still 6 points to interpolate
                        FinalValues = FinalValues(2:end-1, :); %Remove the first and last point
                        
                    elseif length(FinalValues(:,1)) == 5 %If there are 5 points to interpolate
                        switch Case
                            case "X"
                                First = abs(FinalValues(1,2)-X_Cross(j,2));  %Find the distance between the crosspoing and the endpoints of FInalValues
                                Last =  abs(FinalValues(end,2)-X_Cross(j,2));
                            case "Y"
                                First = abs(FinalValues(1,1)-X_Cross(j,1));
                                Last =  abs(FinalValues(end,1)-X_Cross(j,1));
                        end
                        if First>=Last
                            FinalValues = FinalValues(2:end, :); %Remove the first point
                        else
                            FinalValues = FinalValues(1:end-1,:); %Remove the last point
                        end
                    end
                    
                    if length(FinalValues(:,1))>1                       %If we have more than one value to interpolate with 
                        
                        switch Case
                            case "X"
                                Coeff_Ci = obj.LagrangeIP(FinalValues(:,2), X_Cross(j,2));          % Calculate Lagrange IP Coefficents for Ci, Cii, or Ciii 
                            case "Y"
                                Coeff_Ci = obj.LagrangeIP(FinalValues(:,1), X_Cross(j,1));          % Calculate Lagrange IP Coefficents for Ci, Cii, or Ciii 
                        end
%                         if abs(sum(Coeff_Ci)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                               "Coeff_Ci Greater than 1"
%                               %pause
%                          end

                        FinalValues(:,3) = Coeff_Ci*(Coeff_CBar(j));            %Store coefficient values
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;FinalValues]; %Store coordinate with coefficient
                    elseif length(FinalValues(:,1))==1
% %                         "Only One Coefficent"
                        FinalValues(:,3) = Coeff_CBar(j);
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;FinalValues]; %Store coordinate with coefficient
                    else
% %                         "Not enough interpolation points"
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;[0,0,0]];
                    end   
                end           
                
            end
          
% %             if abs(sum(MtrixCorrection_CBar(:,3))+sum(CoeffRho_CBar)-1)>obj.tol
% %               "Coefficient Sum Greater than 1 before fix"
% %              % pause
% %             end
          
          MtrixCorrection_CBar2 = MtrixCorrection_CBar(abs(MtrixCorrection_CBar(:,3))>obj.tol,:); %Remove all points whose coefficents is less than tolerance
          
% %           if abs(sum(MtrixCorrection_CBar2(:,3))+sum(CoeffRho_CBar)-1)>obj.tol
% % % %               "Coefficient Sum Greater than 1"
% %              % pause
% %           end
          out1 = MtrixCorrection_CBar2;
          out2 = CoeffRho_CBar;
        end
        %% Find the location of C_i, C_ii, and C_iii and return the coordinates & coefficients
        function [out1, out2] = Find_CBar_Old(obj, xr_bar, C_Gamma, X_Cross, XSize, poly, CenterID, ka, kd, D, Case)
            %%% NOTE: There is a flaw in this function. If the value C_i,
            %%% C_ii, or C_iii lies on a gridpoint, I do not test for this.
            %%% Instead, I interpolate and if I do not have enough points I
            %%% set to zero. This can be a problem in some cases and should
            %%% be corrected. May 29, 2021.
            
            
            if XSize(1)==0                                                         % if there are no points in X or only one point in X, do nothing
                CoeffRho_CGamma = 0;
                CoeffRho_CBar   = 0;
                MtrixCorrection_CBar = [0,0,0];
                
            else 
                if sqrt(abs(C_Gamma(1,1)-X_Cross(1,1)).^2+abs(C_Gamma(1,2)-X_Cross(1,2)).^2)>sqrt(abs(C_Gamma(1,1)-X_Cross(end,1)).^2+abs(C_Gamma(1,2)-X_Cross(end,2)).^2)
                    X_Cross = flip(X_Cross);
                end
                
                MtrixCorrection_CBar = [];
                xi_r = sqrt(min(abs(C_Gamma(1,1)-X_Cross(:,1))).^2+min(abs(C_Gamma(1,2)-X_Cross(:,2))).^2);       % Find the distance between x_i and x_r -> Fixed on June 4
                if XSize(1)==1
                    xi_ii = 0;
                else 
                    xi_ii   = sqrt(abs(X_Cross(2,1)-X_Cross(1,1)).^2+abs(X_Cross(2,2)-X_Cross(1,2)).^2);                                                                    % Find the distance between x_i and x_r
                end
                
                Coeff_CGamma    = obj.LagrIP_CGamma(XSize(1), xi_r, xi_ii, ka, kd, D);                    % Calculate Lagrange IP coefficents for C_Gamma
                CoeffRho_CGamma = Coeff_CGamma(end);                                              % Identify the coefficent on the rho term
                
%                 if abs(sum(Coeff_CGamma)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                       "Coeff_Ci Greater than 1"
%                       pause
%                 end
                
                if abs(xr_bar)>obj.tol
                    Coeff_CBar = obj.LagrIP_CBar(XSize(1), xr_bar, xi_r, xi_ii, Coeff_CGamma); % Calculate Lagrange IP coefficents for C_bar
                else 
                    Coeff_CBar = Coeff_CGamma;                                                  % If C_bar = C_Gamma, use the same Lagrange IP
                end 
                CoeffRho_CBar = Coeff_CBar(end); %Identify coefficient on rho term
%                 if abs(sum(Coeff_CBar)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                       "Coeff_Ci Greater than 1"
%                      % pause
%                 end
                
                switch Case
                    case "X"
                        Y_Below1=X_Cross(:,2)-mod(X_Cross(:,2),obj.dx);                                          % Find the gridpoint below or to the left of each cross point
                        Y_Vals=[Y_Below1-2*obj.dx,Y_Below1-obj.dx,Y_Below1,Y_Below1+obj.dx,Y_Below1+2*obj.dx,Y_Below1+3*obj.dx]; % Find six points in x surrounding each point in y 
                        X_Vals = Y_Vals*0 + X_Cross(:,1);
                        
                    case "Y"
                        Y_Below1=X_Cross(:,1)-mod(X_Cross(:,1),obj.dx);                                          % Find the gridpoint below or to the left of each cross point
                        X_Vals=[Y_Below1-2*obj.dx,Y_Below1-obj.dx,Y_Below1,Y_Below1+obj.dx,Y_Below1+2*obj.dx,Y_Below1+3*obj.dx]; % Find six points in x surrounding each point in y
                        Y_Vals = X_Vals*0 + X_Cross(:,2);
                end

                for j=1:XSize(1)
                    
                    [Ind, ~] = inpolygon(X_Vals(j,:),Y_Vals(j,:),poly(CenterID).X,poly(CenterID).Y); %Identify which points fall within the boundary
                    
                    FinalValues = [X_Vals(j,:)',Y_Vals(j,:)',zeros(length(Y_Vals(j,:)),1)]; %Store values as a matrix of form [X,Y]
                    FinalValues=FinalValues(Ind'==1 & FinalValues(:,1)>obj.minX & FinalValues(:,1)<obj.maxX & FinalValues(:,2)>obj.minY & FinalValues(:,2)<obj.maxY,:);  %Remove all points that fall outside the domain of interest or on the boundary
                    
                    CheckDist = find(abs(FinalValues(1:end-1)-FinalValues(2:end))-obj.dx)>obj.tol;
                    
                    if ~isempty(CheckDist) & length(FinalValues)>1 %Check to see if a point falls into another domain
                        
                    end
                    
                    
                    if length(FinalValues(:,1)) == 6 %If there are still 6 points to interpolate
                        FinalValues = FinalValues(2:end-1, :); %Remove the first and last point
                        
                    elseif length(FinalValues(:,1)) == 5 %If there are 5 points to interpolate
                        switch Case
                            case "X"
                                First = abs(FinalValues(1,2)-X_Cross(j,2));  %Find the distance between the crosspoing and the endpoints of FInalValues
                                Last =  abs(FinalValues(end,2)-X_Cross(j,2));
                            case "Y"
                                First = abs(FinalValues(1,1)-X_Cross(j,1));
                                Last =  abs(FinalValues(end,1)-X_Cross(j,1));
                        end
                        if First>=Last
                            FinalValues = FinalValues(2:end, :); %Remove the first point
                        else
                            FinalValues = FinalValues(1:end-1,:); %Remove the last point
                        end
                    end
                    
                    if length(FinalValues(:,1))>1                       %If we have more than one value to interpolate with 
                        
                        switch Case
                            case "X"
                                Coeff_Ci = obj.LagrangeIP(FinalValues(:,2), X_Cross(j,2));          % Calculate Lagrange IP Coefficents for Ci, Cii, or Ciii 
                            case "Y"
                                Coeff_Ci = obj.LagrangeIP(FinalValues(:,1), X_Cross(j,1));          % Calculate Lagrange IP Coefficents for Ci, Cii, or Ciii 
                        end
%                         if abs(sum(Coeff_Ci)-1)>obj.tol %Check that the sum of the coefficients is not greater than 1
%                               "Coeff_Ci Greater than 1"
%                               %pause
%                          end

                        FinalValues(:,3) = Coeff_Ci*(Coeff_CBar(j));            %Store coefficient values
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;FinalValues]; %Store coordinate with coefficient
                    elseif length(FinalValues(:,1))==1
% %                         "Only One Coefficent"
                        FinalValues(:,3) = Coeff_CBar(j);
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;FinalValues]; %Store coordinate with coefficient
                    else
% %                         "Not enough interpolation points"
                        MtrixCorrection_CBar=[MtrixCorrection_CBar;[0,0,0]];
                    end   
                end           
                
            end
          
% %             if abs(sum(MtrixCorrection_CBar(:,3))+sum(CoeffRho_CBar)-1)>obj.tol
% %               "Coefficient Sum Greater than 1 before fix"
% %              % pause
% %             end
          
          MtrixCorrection_CBar2 = MtrixCorrection_CBar(abs(MtrixCorrection_CBar(:,3))>obj.tol,:); %Remove all points whose coefficents is less than tolerance
          
% %           if abs(sum(MtrixCorrection_CBar2(:,3))+sum(CoeffRho_CBar)-1)>obj.tol
% % % %               "Coefficient Sum Greater than 1"
% %              % pause
% %           end
          out1 = MtrixCorrection_CBar2;
          out2 = CoeffRho_CBar;
        end
        
        %% Create the polygon for each subdomain and identify which discretization point lies in each domain.
        function [out] = CreatePolygon(obj, poly, Case)
        % X and Y should be a mesh created based on the size of the domain and the
        % spacing dx and dy.
            SizeX=size(obj.Xm);
            SizeY=size(obj.Ym);
            xm=reshape(obj.Xm',SizeX(1)*SizeX(2),1); %Use prime to ensure I reshape by row rather than column
            ym=reshape(obj.Ym',SizeY(1)*SizeY(2),1);
            
            SizeXNB = size(obj.Xm(2:end-1,2:end-1)); %No boundaries
            SizeYNB = size(obj.Ym(2:end-1,2:end-1)); %No boundaries
            xmNB = reshape(obj.Xm(2:end-1,2:end-1)',SizeXNB(1)*SizeXNB(2),1); %Use prime to ensure I reshape by row rather than column
            ymNB = reshape(obj.Ym(2:end-1,2:end-1)',SizeYNB(1)*SizeYNB(2),1);
            
            
            for i = 1:obj.NumDomains                             %Given each set of (X,Y) coordinates that form a row in poly
                poly(i).pgon = polyshape(poly(i).X,poly(i).Y); %Create polygon from vertices
                [in]=inpolygon(xm,ym,poly(i).X,poly(i).Y);     %Check which points are in the domain of interest
                poly(i).Locate = find(in==1);
                poly(i).InX = xm(poly(i).Locate);                 %Save the set of (X,Y) points in each subdomain
                poly(i).InY = ym(poly(i).Locate);
                
                [inNB]=inpolygon(xmNB,ymNB,poly(i).X,poly(i).Y);     %Check which points are in the domain of interest
                poly(i).LocateNB = find(inNB==1);

                poly(i).InXNB = xmNB(poly(i).LocateNB);                 %Save the set of (X,Y) points in each subdomain without the boundries
                poly(i).InYNB = ymNB(poly(i).LocateNB);
                
            end

            switch Case
                case "Plot"
                    figure
                    for i = 1:obj.NumDomains
                        if i == 1
                            plot(poly(1).pgon,'FaceColor','red')
                            hold on
                        elseif i == 2
                            plot(poly(2).pgon,'FaceColor','blue')
                            hold on
                        elseif i == 3
                            plot(poly(3).pgon,'FaceColor','green')
                            hold on
                        else
                            "More than three subdomains"
                            plot(obj.poly(3).pgon,'FaceColor','cyan')
                            hold on
                        end
                    end
                    xlabel("X")
                    ylabel("Y")
                    title("Plot of Subdomains. Red is Domain 1, Blue is Domain 2, Green is Domain 3, and all others are Cyan.")

                case "No Plot"
                    "No plot of subdomains requested"
            end
            out = poly;
        end
        
        %% Create the network for each subdomain and identify which discretization points line in each channel.
        function [out] = CreateMatrixNetwork(obj, k, Case)
            
            if obj.NumChannels == 1
               switch Case
                   case "ZeroDirichlet"
                       i = obj.NetMtrix(1).Size; %Set the size of the matrix based on the number of points in the discretization minus two for the BC
                       r = obj.NetMtrix(1).r*k; %Update the coefficient to include time
                       
                       off(1:i-3)=r; %Create the off diagonals of the matrix
                       onL(1:i-2)=1+2*r; %Create the main diagonals of the right and left matrix
                       onR(1:i-2)=1-2*r;

                       Mtrix.LHS=diag(onL,0)+diag(-off,-1)+diag(-off,1); % Create the centered spatial difference matrix
                       Mtrix.RLHS=diag(onR,0)+diag(off,-1)+diag(off,1); 
                       Mtrix.RhoID_X = obj.NetMtrix(1).DiscIn2D_X(2:end-1,1);
                       Mtrix.RhoID_Y = obj.NetMtrix(1).DiscIn2D_Y(2:end-1,1);
                       Mtrix.ID = ones(obj.NetMtrix(1).Size,1);
      
               end
               
            elseif obj.NumChannels == 2
                switch Case
                   case "ZeroDirichlet"
                       
                       Size = [obj.NetMtrix(1).Size, obj.NetMtrix(2).Size]; %Find the number of discretization points in channels 1 and 2

                       r = [obj.NetMtrix(1).r*k, obj.NetMtrix(2).r*k]; %Correct the coefficients for channels 1 and 2
                       
                      LocVertex1 = obj.NetMtrix(1).LocVertex;
                      LocVertex2 = obj.NetMtrix(2).LocVertex;
                       
                       Fix = sparse(1,Size(1)-2+Size(2)-2);
                       Coeff = 1/((((obj.NetMtrix(1).D/obj.NetMtrix(1).drho) + (obj.NetMtrix(2).D)/obj.NetMtrix(2).drho))*(3/2)); %Find the correct coefficient
                       
                       if (LocVertex2==1 & LocVertex1==Size(1)) %If channel 2 is connected on the right to channel 1

                           Fix(Size(1)-2) =  Coeff*2*obj.NetMtrix(1).D/obj.NetMtrix(1).drho; %Create the row of coefficients needed to edit
                           Fix(Size(1)-3) = -Coeff*(1/2)*obj.NetMtrix(1).D/obj.NetMtrix(1).drho;
                           Fix(Size(1)-1) =  Coeff*2*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           Fix(Size(1))   = -Coeff*(1/2)*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           
                           FixLoc = [Size(1)-2;Size(1)-1];
                           
                       elseif (LocVertex1==1 & LocVertex2==Size(2)) %If channel 2 is connected on the right to channel 1
                           
                           Fix(1) =  Coeff*2*obj.NetMtrix(1).D/obj.NetMtrix(1).drho; %Create the row of coefficients needed to edit
                           Fix(2) = -Coeff*(1/2)*obj.NetMtrix(1).D/obj.NetMtrix(1).drho;
                           Fix(end) =  Coeff*2*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           Fix(end-1)   = -Coeff*(1/2)*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           
                           FixLoc = [1;Size(1)-2+Size(2)-2];
                           
                       elseif (LocVertex2==1 & LocVertex1==1) %If channel 2 is connected on the right to channel 1
                           
                           Fix(1) =  Coeff*2*obj.NetMtrix(1).D/obj.NetMtrix(1).drho; %Create the row of coefficients needed to edit
                           Fix(2) = -Coeff*(1/2)*obj.NetMtrix(1).D/obj.NetMtrix(1).drho;
                           Fix(Size(1)-1) =  Coeff*2*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           Fix(Size(1))   = -Coeff*(1/2)*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           
                           FixLoc = [1; Size(1)-1];
                           
                       elseif (LocVertex2==Size(2) & LocVertex1==Size(1)) %If channel 2 is connected on the right to channel 1
                           
                           Fix(Size(1)-2) =  Coeff*2*obj.NetMtrix(1).D/obj.NetMtrix(1).drho; %Create the row of coefficients needed to edit
                           Fix(Size(1)-3) = -Coeff*(1/2)*obj.NetMtrix(1).D/obj.NetMtrix(1).drho;
                           Fix(end) =  Coeff*2*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           Fix(end-1)   = -Coeff*(1/2)*obj.NetMtrix(2).D/obj.NetMtrix(2).drho;
                           
                           FixLoc = [Size(1)-1;Size(1)-2+Size(2)-2];
                       end
                       
                       i = Size(1); %For the first channel
                       off(1:i-3)=r(1); %Create the off diagonals of the matrix
                       onL(1:i-2)=1+2*r(1); %Create the main diagonals of the right and left matrix
                       onR(1:i-2)=1-2*r(1);

                       LHS_1=diag(onL,0)+diag(-off,-1)+diag(-off,1); % Create the centered spatial difference matrix
                       RLHS_1=diag(onR,0)+diag(off,-1)+diag(off,1); 

                       i = Size(2); %For the second channel
                       off(1:i-3)=r(2); %Create the off diagonals of the matrix
                       onL(1:i-2)=1+2*r(2); %Create the main diagonals of the right and left matrix
                       onR(1:i-2)=1-2*r(2);

                       LHS_2=diag(onL,0)+diag(-off,-1)+diag(-off,1); % Create the centered spatial difference matrix
                       RLHS_2=diag(onR,0)+diag(off,-1)+diag(off,1); 

                       LHS = blkdiag(LHS_1,LHS_2); %Create full matrix
                       RLHS = blkdiag(RLHS_1,RLHS_2);

                       LHS(FixLoc,:)  = LHS(FixLoc,:) -[r(1)*Fix;r(2)*Fix]; %Update matrix with fix times the correct r for the domain
                       RLHS(FixLoc,:) = RLHS(FixLoc,:)+[r(1)*Fix;r(2)*Fix];

                       Mtrix.RLHS(:,:) = RLHS;
                       Mtrix.LHS(:,:) = LHS;
                       Mtrix.Junction(:,1) = Fix;
                       Mtrix.RhoID_X = [obj.NetMtrix(1).DiscIn2D_X(2:end-1,1);obj.NetMtrix(2).DiscIn2D_X(2:end-1,1)];
                       Mtrix.RhoID_Y = [obj.NetMtrix(1).DiscIn2D_Y(2:end-1,1);obj.NetMtrix(2).DiscIn2D_Y(2:end-1,1)];
                       Mtrix.ID = [ones(Size(1)-2,1);ones(Size(2)-2,1)+1];
                end
               
                elseif obj.NumChannels == 3
                    
                    "Not Written Yet"
            end
            
         out = Mtrix;    
       
        end
        
        %% Create the matrix of coefficents for each subdomain
        function [out] = CreateMatrixFields(obj, poly, dt,Case, Case2)
            
            for i = 1:obj.NumDomains %For each subdomain
                Index = poly(i).Locate; %Set the location of the points in domain i
                InX = poly(i).InXNB;      %Set the coordinates of the points in domain i
                InY = poly(i).InYNB;
                [InXOnBound, InYOnBound] = CheckBoundary(obj,InX,InY); %Find points on the outer boundary of computational domain
                R = (poly(i).D*dt)/(2*(obj.dx^2));
                switch Case  
                    case "ZeroDirichlet"
                        PntsOnBoundary = unique([InXOnBound; InYOnBound]);  %Find the uniquely defined points on the boundary
                        InX(PntsOnBoundary)=[];                             %Remove points on the boundary from the points of interest
                        InY(PntsOnBoundary)=[];
                        Index(PntsOnBoundary)=[];
                        Len = length(InX);
                        RLHS = sparse(Len,Len);
                        Laplace = sparse(Len,Len);
                        InterpRho = [];  %Create holding matrix for C_Gamma and the coefficient
                        LHS = sparse(Len,Len);%Create base matrix of coefficents 
                        
                        for j = 1:Len                                               %For each remaining point in the domain:
                            LocCGammaX = [];
                            LocCGammaY = [];
                            NeigborsX = [];
                            NeigborsY = []; 
                            
                                try
                                    Loc1 = find(abs(obj.Xm(1,:) - (InY(j)-obj.dy))<obj.tol); 
                                    if ~isempty(Loc1)
                                        NeigborsY = [obj.Xm(1,Loc1(1))];
                                        NeigborsX = [InX(j)]; 
                                    end
                                catch
                                end
                                
                                try
                                    Loc2 = find(abs(obj.Xm(1,:) - (InY(j)+obj.dy))<obj.tol);
                                    if ~isempty(Loc2)
                                        NeigborsY = [NeigborsY,obj.Xm(1,Loc2(1))];
                                        NeigborsX = [NeigborsX,InX(j)]; 
                                    end
                                catch
                                end
                                
                                try
                                    Loc2 = find(abs(obj.Xm(1,:) - (InX(j)+obj.dx))<obj.tol);
                                    if ~isempty(Loc2)
                                        NeigborsX = [NeigborsX,obj.Xm(1,Loc2(1))];
                                        NeigborsY = [NeigborsY,InY(j)];
                                    end
                                catch     
                                end
                                
                                try
                                    Loc1 = find(abs(obj.Xm(1,:) - (InX(j)-obj.dx))<obj.tol); 
                                    if ~isempty(Loc1)
                                        NeigborsX = [NeigborsX,obj.Xm(1,Loc1(1))];
                                        NeigborsY = [NeigborsY,InY(j)];
                                    end
                                catch     
                                end
                                
%                                 NeigborsX = [InX(j),InX(j),InX(j+1),InX(j-1)];            %Find the four neighbors [Above,Below,Right,Left]
%                                 Loc1 = find(abs(InY - (InY(j)-obj.dy))<obj.tol);              %Find the location of the neighbor in InY. If I just add or subtract dy, the inpolygon function may not work
%                                 Loc2 = find(abs(InY - (InY(j)+obj.dy))<obj.tol);              %Find the location of the neighbor in InY. If I just add or subtract dy, the inpolygon function may not work
%                                 NeigborsY = [InY(Loc2(1)),InY(Loc1(1)),InY(j),InY(j)];

                          
                            
                            
                            
                            [OutXBound, OutYBound] = CheckBoundary(obj,NeigborsX,NeigborsY);    %Find points on the outer boundary of computational domain
                            PntsOnBoundary = unique([OutXBound; OutYBound]);                    %Find the uniquely defined points on the boundary
                            NeigborsX(PntsOnBoundary)=[];                                       %Remove points on the boundary from the points of interest
                            NeigborsY(PntsOnBoundary)=[];

                            [InP,~] =  inpolygon(NeigborsX,NeigborsY,poly(i).X,poly(i).Y);     %Check for neighbors that fall into another domain
                            LocInOtherDomain = find(InP==0);
                            if ~isempty(LocInOtherDomain)                                       %If some neighbors fall into another domain, estimate their values
%                                 plot(InX(j),InY(j), "g*")
%                                 hold on 
                                NeigborXOtherDomain = NeigborsX(LocInOtherDomain);                 %Separate out neighbors in domain vs those out of the domain
                                NeigborYOtherDomain = NeigborsY(LocInOtherDomain);
                                NeigborsX(LocInOtherDomain) = [];
                                NeigborsY(LocInOtherDomain) = [];

                                
                                for k = 1:length(NeigborXOtherDomain)                              %For each neighbor crossing into another domain
                                    [EBMCoeff, Coeff_Rho, C_Gamma, C_i, LineID] = obj.EBM(poly, i, [InX(j),InY(j)],[NeigborXOtherDomain(k), NeigborYOtherDomain(k)]); %Use the Embedded Boundary Method to estimate the concentration at the crosspoint
                                    switch Case2
                                        case "Plot" %Plot results of the EBM
                                            obj.PlotEBM(poly, i, [InX(j),InY(j)],[NeigborXOtherDomain(k), NeigborYOtherDomain(k)], C_Gamma, C_i, EBMCoeff)
                                    end 
                                    
                                    
                                    if ~isempty(EBMCoeff) %Check that EBM is not empty
                                        Size = size(EBMCoeff); %Find the number of points used to interpolate C_Bar
                                        InterpRho = [InterpRho;j,C_Gamma(1), C_Gamma(2), R*Coeff_Rho, LineID]; %Stores the [location in matrix, Location of C_Gamma in X and Y, Coefficient for that point, Line ID]
                                        for l = 1:Size(1)
                                            Loc = find(abs(EBMCoeff(l,1)-InX)<obj.tol &  abs(EBMCoeff(l,2)-InY)<obj.tol);   %Find the location of each neighbor in the vector of all points
                                           % RLHS(j,Loc) = RLHS(j,Loc) + R*EBMCoeff(l,3);                                                        %Add the correct coefficent in that spot of the matrix
                                           % LHS(j,Loc) = LHS(j,Loc) + R*EBMCoeff(l,3);
                                            Laplace(j,Loc) = Laplace(j,Loc) +R*EBMCoeff(l,3);
                                        end
                                        
                                        "A point falls in another domain and requires the EBM"; 
                                    else
                                        "Empty";
                                    end
                                end
                            end

                            for k = 1:length(NeigborsX)                                                         %For each of the remaining neighbors in the domain of interest
                                Loc = find(abs(NeigborsX(k)-InX)<obj.tol &  abs(NeigborsY(k)-InY)<obj.tol);         %Find the location of each neighbor in the vector of all points
                               % RLHS(j,Loc) = RLHS(j,Loc) + R;                                                        %Add the correct coefficent in that spot of the matrix
                               % LHS(j,Loc) = LHS(j,Loc) + R;
                                Laplace(j,Loc) = Laplace(j,Loc) +R;
                            end
                            Laplace(j,j) = Laplace(j,j) - 4*R;
                           % RLHS(j,j) = 1+(RLHS(j,j) - 4*R);                                          %Add the correct coefficent at the center of the stencil
                           % LHS(j,j) = 1-(LHS(j,j) - 4*R); 
                        end           
                end
                Diag = ones(Len,1);
                Matrix1(i).Laplace = Laplace;
                Matrix1(i).RLHS = spdiags(Diag,0,Len,Len)+Laplace;
                Matrix1(i).LHS = spdiags(Diag,0,Len,Len)-Laplace;   
                Matrix1(i).CGamma = InterpRho;
            end
            out = Matrix1; 
        end
        
        %% Find concentration at C_bar using Embedded Boundary Method
        function [out1, out2, out3, out4, out5] = EBM(obj,poly, CenterID, Center, Cross)
            
            LocCenter = find(obj.Network.Boundary(:,1)==CenterID); %Find the domain in which the center lives
            
            NeigborID = unique(obj.Network.Boundary(LocCenter,2)); %Find the ID of the neighbors
            
            for i = 1:length(NeigborID)
                [InP,~] =  inpolygon(Cross(1),Cross(2),poly(NeigborID(i)).X,poly(NeigborID(i)).Y);     %Check which domain the Crosspoint falls in
                if InP == 1
                    CrossID = NeigborID(i); %Store the domain ID
                else
                    "Error with Initial Network Setup"
                    pause
                end
            end
            
            BndryLoc = find(obj.Network.Boundary(:,1)==CenterID & obj.Network.Boundary(:,2)==CrossID | obj.Network.Boundary(:,2)==CenterID & obj.Network.Boundary(:,1)==CrossID);  %Find the correct row(s) where the boundary information is stored
            LineID = 0;
            xi_bar = 100000;
            for i = 1:length(BndryLoc)
                Start = [obj.Network.vertexX(BndryLoc(i),1),obj.Network.vertexY(BndryLoc(i),1)]; %Identify the (x,y) pair where the boundary starts
                End   = [obj.Network.vertexX(BndryLoc(i),2),obj.Network.vertexY(BndryLoc(i),2)]; %Identify the (x,y) pair where the boundary ends
                Check = obj.FindChannel(Start, End, Center, Cross); %Check if the line falls between the crosspoint and center of stencil
               
                if Check == 1
                    [xi_bar_Temp, C_Gamma_Temp, Type_Temp, m_and_b_Temp] = obj.point_to_line(Cross, Start, End, "No Plot"); %Find the normal passing through the crosspoint and normal to a boundary line
                    if xi_bar_Temp < xi_bar                                        %Find if the distance from C_Bar to the line is less than the current option
                        xi_bar = xi_bar_Temp;                                      %If yes, x_bar_Temp becomes the new minimum 
                        C_Gamma = C_Gamma_Temp;
                        Type = Type_Temp; 
                        m_and_b = m_and_b_Temp;
                        ka = obj.Network.ka(i);
                        kd = obj.Network.kd(i);
                        LineID = obj.Network.LineID(BndryLoc(i));
                    end 
                end
            end
            
            [X_i, Y_i, CaseGamma] = obj.GenerateStencil(Center, Cross, C_Gamma, Type, m_and_b(1)); %Find the neigboring gridlines in X and Y and whether the flux on the BC should be negative or positive
            if length(X_i)>6 
                'X_i is longer than 6!'
                pause
            end
            switch Type
                case "Horizontal"
                    X_CrossT=[0,0;0,0;0,0;0,0;0,0;0,0]; %Find C_i, C_ii, C_iii 
                    Y_CrossT=[Center(1)*ones(length(Y_i),1), Y_i'];
                case "Angled"
                    X_CrossT=[X_i',obj.NormalLine(m_and_b(1), m_and_b(2),X_i,"y")']; %Find C_i, C_ii, C_iii 
                    Y_CrossT=[obj.NormalLine(m_and_b(1), m_and_b(2),Y_i,"x")', Y_i'];
            end
            
            RightId = 1;
            LeftId = 1;
            RightIndex = [3,2,1];
            LeftIndex = [4,5,6];
            X_Cross = X_CrossT;
            Index = 1;
            Count = 0;
            
            while Index < 4 
                if RightId ==1
                    [InX_Cross,~] = inpolygon(X_CrossT(RightIndex(Index),1),X_CrossT(RightIndex(Index),2),poly(CenterID).X,poly(CenterID).Y); %Identify which points (C_i, C_ii, C_iii) live in the domain of interest
                  
                    if InX_Cross==1 & abs(X_CrossT(RightIndex(Index),1)-obj.minX)>obj.tol & abs(X_CrossT(RightIndex(Index),1)-obj.maxX)>obj.tol & abs(X_CrossT(RightIndex(Index),2)-obj.minY)>obj.tol & abs(X_CrossT(RightIndex(Index),2)-obj.maxY)>obj.tol
                        Count = Count+1;
                        if Count > 3
                            RightId = 0;
                            X_CrossT(RightIndex(Index),:)=1000000;
                        end       
                    else
                        RightId = 0;
                        X_CrossT(RightIndex(Index),:)=1000000;
                    end
                else
                    X_CrossT(RightIndex(Index),:)=1000000;
                end
                
                if LeftId ==1
                    [InX_Cross,~] = inpolygon(X_CrossT(LeftIndex(Index),1),X_CrossT(LeftIndex(Index),2),poly(CenterID).X,poly(CenterID).Y); %Identify which points (C_i, C_ii, C_iii) live in the domain of interest
                    if InX_Cross==1 & abs(X_CrossT(LeftIndex(Index),1)-obj.minX)>obj.tol & abs(X_CrossT(LeftIndex(Index),1)-obj.maxX)>obj.tol & abs(X_CrossT(LeftIndex(Index),2)-obj.minY)>obj.tol & abs(X_CrossT(LeftIndex(Index),2)-obj.maxY)>obj.tol
                        Count = Count + 1;
                        if Count > 3
                            LeftId = 0;
                            X_CrossT(LeftIndex(Index),:)=1000000;
                        end 
                    else
                        LeftId = 0;
                        X_CrossT(LeftIndex(Index),:)=1000000;
                    end
                else
                    X_CrossT(LeftIndex(Index),:)=1000000;
                end
                
                Index = Index+1;
                
            end
            
            RightId = 1;
            LeftId = 1;
            RightIndex = [3,2,1];
            LeftIndex = [4,5,6];
            Y_Cross = Y_CrossT;
            Index = 1;
            Count = 0;
            
            while Index < 4
                if RightId ==1
                   [InY_Cross,~] = inpolygon(Y_CrossT(RightIndex(Index),1),Y_CrossT(RightIndex(Index),2),poly(CenterID).X,poly(CenterID).Y);
                   if InY_Cross==1 & abs(Y_CrossT(RightIndex(Index),1)-obj.minX)>obj.tol & abs(Y_CrossT(RightIndex(Index),1)-obj.maxY)>obj.tol & abs(Y_CrossT(RightIndex(Index),2)-obj.minY)>obj.tol & abs(Y_CrossT(RightIndex(Index),2)-obj.maxY)>obj.tol
                        Count = Count + 1;
                        if Count > 3
                            RightId = 0;
                            Y_CrossT(RightIndex(Index),:)=1000000;
                            'Max points reached'
                        end 
                    else
                        RightId = 0;
                        Y_CrossT(RightIndex(Index),:)=1000000;
                    end
                else
                   Y_CrossT(RightIndex(Index),:)=1000000;
                end
                
                if LeftId ==1
                  [InY_Cross,~] = inpolygon(Y_CrossT(LeftIndex(Index),1),Y_CrossT(LeftIndex(Index),2),poly(CenterID).X,poly(CenterID).Y);
                    if InY_Cross==1 & abs(Y_CrossT(LeftIndex(Index),1)-obj.minX)>obj.tol & abs(Y_CrossT(LeftIndex(Index),1)-obj.maxY)>obj.tol & abs(Y_CrossT(LeftIndex(Index),2)-obj.minY)>obj.tol & abs(Y_CrossT(LeftIndex(Index),2)-obj.maxY)>obj.tol
                        Count = Count + 1;
                        if Count > 3
                            LeftId = 0;
                            Y_CrossT(LeftIndex(Index),:)=1000000;
                        end 
                    else
                        LeftId = 0;
                        Y_CrossT(LeftIndex(Index),:)=1000000;
                    end
                else
                    Y_CrossT(LeftIndex(Index),:)=1000000;
                end
                Index = Index+1;
                
            end
            X_Cross = X_Cross(X_CrossT(:,1)<1000000 & X_CrossT(:,2)<1000000,:);
            Y_Cross = Y_Cross(Y_CrossT(:,1)<1000000 & Y_CrossT(:,2)<1000000,:);
            XSize=size(X_Cross); % Find the number of points that cross in X and Y
            YSize=size(Y_Cross);
            
            
            if isempty(X_Cross)
                dX1 = 1000000;
                dX2 = 1000000;
            else
                dX1 = sqrt((C_Gamma(1)-X_Cross(1,1))^2+(C_Gamma(2)-X_Cross(1,2))^2);
                dX2 = sqrt((C_Gamma(1)-X_Cross(end,1))^2+(C_Gamma(2)-X_Cross(end,2))^2);
            end
            
            
            if isempty(Y_Cross)
                dY1 = 1000000;
                dY2 = 1000000;
            else   
                dY1 = sqrt((C_Gamma(1)-Y_Cross(1,1))^2+(C_Gamma(2)-Y_Cross(1,2))^2);
                dY2 = sqrt((C_Gamma(1)-Y_Cross(end,1))^2+(C_Gamma(2)-Y_Cross(end,2))^2);
            end
            
            switch CaseGamma
                case "Equal"
                    switch Type
                        case "Horizontal"
                            if Center(2)>Y_Cross(1,2)
                                Y_Cross = [Center;Y_Cross];
                            else
                                Y_Cross = [Y_Cross;Center];
                            end
                            InterpCBar = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
                            Coeff_Rho = 0;
                            out4 = Y_Cross;

%                             if Center(2)>Y_Cross(1,2)
%                                 Y_Cross = [Center;Y_Cross];
%                                 InterpCBar_Temp = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
%                                 Coeff_Rho = InterpCBar_Temp(1,3);
%                                 InterpCBar = InterpCBar_Temp(2:end,:);
%                                 out4 = Y_Cross(2:end,:);
%                             else
%                                 Y_Cross = [Y_Cross;Center];
%                                 InterpCBar_Temp = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
%                                 Coeff_Rho = InterpCBar_Temp(end,3);
%                                 InterpCBar = InterpCBar_Temp(1:end-1,:);
%                                 out4 = Y_Cross(1:end-1,:);
%                             end
                            
                            
%                               if Center(2)>Cross(2)
%                                   InterpCBar = [Y_Cross,[1;(2*obj.dx*ka)/poly(CenterID).D]];
%                                   Coeff_Rho = -(2*obj.dx*kd)/poly(CenterID).D;
%                                   out4 = Y_Cross;
%                               else
%                                   InterpCBar = [Y_Cross,[-(2*obj.dx*ka)/poly(CenterID).D;1]];
%                                   Coeff_Rho = (2*obj.dx*kd)/poly(CenterID).D;
%                                   out4 = Y_Cross;
%                               end                              

                        case "Angled"
                            if XSize(1)~=0 || YSize(1)~=0                  %If at least one cross point exists in X or Y

                                if max(dX1,dX2)<=max(dY1,dY2) %If there are more points in X OR if there are equal points and the slope of the boundary is greater than or equal to 1
                                   InterpCBar = [X_Cross,obj.LagrangeIP(X_Cross(:,1), Cross(1))'];
                                   Coeff_Rho = 0;
                                   out4 = X_Cross;
                                   
                                elseif max(dX1,dX2)>max(dY1,dY2) %If there are more points in Y OR if there are equal points and the slope of the normal is less than or equal to 1
                                   InterpCBar = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
                                   Coeff_Rho = 0;
                                   out4 = Y_Cross;
                                   
                                else
                                    "No case exists for nonempty X_Cross and Y_Cross"
                                    InterpCBar = [];
                                    Coeff_Rho = 0;
                                    out4 = [0,0];
                                end
                            else 
                                "C_i, C_ii, and C_iii do not lie in the domain";
                                InterpCBar = [];
                                Coeff_Rho = 0;
                                out4 = [0,0];
                            end
                   end
                case "NotEqual"
                    if XSize(1)~=0 || YSize(1)~=0                                  %If at least one cross point exists in X or Y
                        D = poly(CenterID).D;
                        if max(dX1,dX2)<=max(dY1,dY2) %If there are more points in X OR if there are equal points and the slope of the boundary is greater than or equal to 1
%                             
                            [InterpCBar, Coeff_Rho] = obj.Find_CBar(xi_bar, C_Gamma, X_Cross, XSize, poly, CenterID, ka, kd, poly(CenterID).D, "X"); 
                           out4 = X_Cross;
                        elseif max(dX1,dX2)>max(dY1,dY2) %If there are more points in Y OR if there are equal points and the slope of the normal is less than or equal to 1
%                             
                            [InterpCBar, Coeff_Rho] = obj.Find_CBar(xi_bar, C_Gamma, Y_Cross, YSize, poly, CenterID, ka, kd, poly(CenterID).D, "Y");
                           out4 = Y_Cross;
                        else
                            "No case exists for nonempty X_Cross and Y_Cross"
                            InterpCBar = [];
                            Coeff_Rho = 0;
                            out4 = [0,0];
                        end
                    else 
                        "C_i, C_ii, and C_iii do not lie in the domain";
                        InterpCBar = [];
                        Coeff_Rho = 0;
                        out4 = [0,0];
                    end
            end
            out1 = InterpCBar; 
            out2 = Coeff_Rho;
            out3 = C_Gamma;
            out5 = LineID;
        end
        
        %% Find concentration at C_bar using Embedded Boundary Method
        function [out1, out2, out3, out4, out5] = EBM_Old(obj,poly, CenterID, Center, Cross)
            
            LocCenter = find(obj.Network.Boundary(:,1)==CenterID); %Find the domain in which the center lives
            
            NeigborID = unique(obj.Network.Boundary(LocCenter,2)); %Find the ID of the neighbors
            
            for i = 1:length(NeigborID)
                [InP,~] =  inpolygon(Cross(1),Cross(2),poly(NeigborID(i)).X,poly(NeigborID(i)).Y);     %Check which domain the Crosspoint falls in
                if InP == 1
                    CrossID = NeigborID(i); %Store the domain ID
                else
                    "Error with Initial Network Setup"
                    pause
                end
            end
            
            BndryLoc = find(obj.Network.Boundary(:,1)==CenterID & obj.Network.Boundary(:,2)==CrossID | obj.Network.Boundary(:,2)==CenterID & obj.Network.Boundary(:,1)==CrossID);  %Find the correct row(s) where the boundary information is stored
            LineID = 0;
            xi_bar = 100000;
            for i = 1:length(BndryLoc)
                Start = [obj.Network.vertexX(BndryLoc(i),1),obj.Network.vertexY(BndryLoc(i),1)]; %Identify the (x,y) pair where the boundary starts
                End   = [obj.Network.vertexX(BndryLoc(i),2),obj.Network.vertexY(BndryLoc(i),2)]; %Identify the (x,y) pair where the boundary ends
                Check = obj.FindChannel(Start, End, Center, Cross); %Check if the line falls between the crosspoint and center of stencil
               
                if Check == 1
                    [xi_bar_Temp, C_Gamma_Temp, Type_Temp, m_and_b_Temp] = obj.point_to_line(Cross, Start, End, "No Plot"); %Find the normal passing through the crosspoint and normal to a boundary line
                    if xi_bar_Temp < xi_bar                                        %Find if the distance from C_Bar to the line is less than the current option
                        xi_bar = xi_bar_Temp;                                      %If yes, x_bar_Temp becomes the new minimum 
                        C_Gamma = C_Gamma_Temp;
                        Type = Type_Temp; 
                        m_and_b = m_and_b_Temp;
                        ka = obj.Network.ka(i);
                        kd = obj.Network.kd(i);
                        LineID = obj.Network.LineID(BndryLoc(i));
                    end 
                end
            end
            
            [X_i, Y_i, CaseGamma] = obj.GenerateStencil(Center, Cross, C_Gamma, Type, m_and_b(1)); %Find the neigboring gridlines in X and Y and whether the flux on the BC should be negative or positive
            
            switch Type
                case "Horizontal"
                    X_Cross=[0,0]; %Find C_i, C_ii, C_iii 
                    Y_Cross=[Center(1)*ones(length(Y_i),1), Y_i'];
                case "Angled"
                    X_Cross=[X_i',obj.NormalLine(m_and_b(1), m_and_b(2),X_i,"y")']; %Find C_i, C_ii, C_iii 
                    Y_Cross=[obj.NormalLine(m_and_b(1), m_and_b(2),Y_i,"x")', Y_i'];
            end
            
            
            [InX_Cross,~] = inpolygon(X_Cross(:,1),X_Cross(:,2),poly(CenterID).X,poly(CenterID).Y); %Identify which points (C_i, C_ii, C_iii) live in the domain of interest
            [InY_Cross,~] = inpolygon(Y_Cross(:,1),Y_Cross(:,2),poly(CenterID).X,poly(CenterID).Y);
            
            X_Cross=X_Cross(InX_Cross==1 & abs(X_Cross(:,1)-obj.minX)>obj.tol & abs(X_Cross(:,1)-obj.maxX)>obj.tol & abs(X_Cross(:,2)-obj.minY)>obj.tol & abs(X_Cross(:,2)-obj.maxY)>obj.tol,:);  %Remove all points that fall outside the domain of interest or on the boundary
            Y_Cross=Y_Cross(InY_Cross==1 & abs(Y_Cross(:,1)-obj.minX)>obj.tol & abs(Y_Cross(:,1)-obj.maxY)>obj.tol & abs(Y_Cross(:,2)-obj.minY)>obj.tol & abs(Y_Cross(:,2)-obj.maxY)>obj.tol,:);
            
            
            XSize=size(X_Cross); % Find the number of points that cross in X and Y
            YSize=size(Y_Cross);
            
            if (XSize(1) > 3) 
                Diff = X_Cross(1:end-1,1) - X_Cross(2:end,1);
                LocMax = find(abs(Diff-max(Diff))<obj.tol);
                LocMin = find(abs(Diff-min(Diff))<obj.tol);
                Location = find(abs(Diff-min(Diff))>obj.tol);
                if sum(Diff(1:Location)) > sum(Diff(Location:end))
                    X_Cross = X_Cross(1:Location-1,:);
                elseif sum(Diff(1:Location)) < sum(Diff(Location:end))
                    X_Cross = X_Cross(Location+1:end,:);
                else
                    "XSize Incorrect in EBM"
                pause
                end
                XSize = size(X_Cross);
            end
            
            if (YSize(1) > 3) 
                Diff = Y_Cross(1:end-1,2) - Y_Cross(2:end,2);
                LocMax = find(abs(Diff-max(Diff))<obj.tol);
                LocMin = find(abs(Diff-min(Diff))<obj.tol);
                Location = find(abs(Diff-min(Diff))>obj.tol);
                if sum(Diff(1:Location)) > sum(Diff(Location:end))
                    Y_Cross = Y_Cross(1:Location-1,:);
                elseif sum(Diff(1:Location)) < sum(Diff(Location:end)) %#ok<*AND2>
                    Y_Cross = Y_Cross(Location+1:end,:);
                else
                    if length(Y_Cross(:,1))==6
                        Y_Cross = Y_Cross(1:3,:)
                    end
                   % "XSize Incorrect" %#ok<*NOPRT>
                    %pause
                %pause
                end
                YSize = size(Y_Cross);
            end

            if isempty(X_Cross)
                dX1 = 1000000;
                dX2 = 1000000;
            else
                dX1 = sqrt((C_Gamma(1)-X_Cross(1,1))^2+(C_Gamma(2)-X_Cross(1,2))^2);
                dX2 = sqrt((C_Gamma(1)-X_Cross(end,1))^2+(C_Gamma(2)-X_Cross(end,2))^2);
            end
            
            
            if isempty(Y_Cross)
                dY1 = 1000000;
                dY2 = 1000000;
            else   
                dY1 = sqrt((C_Gamma(1)-Y_Cross(1,1))^2+(C_Gamma(2)-Y_Cross(1,2))^2);
                dY2 = sqrt((C_Gamma(1)-Y_Cross(end,1))^2+(C_Gamma(2)-Y_Cross(end,2))^2);
            end
            
            switch CaseGamma
                case "Equal"
                    switch Type
                        case "Horizontal"
                            if Center(2)>Y_Cross(1,2)
                                Y_Cross = [Center;Y_Cross];
                            else
                                Y_Cross = [Y_Cross;Center];
                            end
                            InterpCBar = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
                            Coeff_Rho = 0;
                            out4 = Y_Cross;

%                             if Center(2)>Y_Cross(1,2)
%                                 Y_Cross = [Center;Y_Cross];
%                                 InterpCBar_Temp = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
%                                 Coeff_Rho = InterpCBar_Temp(1,3);
%                                 InterpCBar = InterpCBar_Temp(2:end,:);
%                                 out4 = Y_Cross(2:end,:);
%                             else
%                                 Y_Cross = [Y_Cross;Center];
%                                 InterpCBar_Temp = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
%                                 Coeff_Rho = InterpCBar_Temp(end,3);
%                                 InterpCBar = InterpCBar_Temp(1:end-1,:);
%                                 out4 = Y_Cross(1:end-1,:);
%                             end
                            
                            
%                               if Center(2)>Cross(2)
%                                   InterpCBar = [Y_Cross,[1;(2*obj.dx*ka)/poly(CenterID).D]];
%                                   Coeff_Rho = -(2*obj.dx*kd)/poly(CenterID).D;
%                                   out4 = Y_Cross;
%                               else
%                                   InterpCBar = [Y_Cross,[-(2*obj.dx*ka)/poly(CenterID).D;1]];
%                                   Coeff_Rho = (2*obj.dx*kd)/poly(CenterID).D;
%                                   out4 = Y_Cross;
%                               end                              

                        case "Angled"
                            if XSize(1)~=0 || YSize(1)~=0                  %If at least one cross point exists in X or Y

                                if max(dX1,dX2)<=max(dY1,dY2) %If there are more points in X OR if there are equal points and the slope of the boundary is greater than or equal to 1
                                   InterpCBar = [X_Cross,obj.LagrangeIP(X_Cross(:,1), Cross(1))'];
                                   Coeff_Rho = 0;
                                   out4 = X_Cross;
                                   
                                elseif max(dX1,dX2)>max(dY1,dY2) %If there are more points in Y OR if there are equal points and the slope of the normal is less than or equal to 1
                                   InterpCBar = [Y_Cross,obj.LagrangeIP(Y_Cross(:,2), Cross(2))'];
                                   Coeff_Rho = 0;
                                   out4 = Y_Cross;
                                   
                                else
                                    "No case exists for nonempty X_Cross and Y_Cross"
                                    InterpCBar = [];
                                    Coeff_Rho = 0;
                                    out4 = [0,0];
                                end
                            else 
                                "C_i, C_ii, and C_iii do not lie in the domain";
                                InterpCBar = [];
                                Coeff_Rho = 0;
                                out4 = [0,0];
                            end
                   end
                case "NotEqual"
                    if XSize(1)~=0 || YSize(1)~=0                                  %If at least one cross point exists in X or Y
                        D = poly(CenterID).D;
                        if max(dX1,dX2)<=max(dY1,dY2) %If there are more points in X OR if there are equal points and the slope of the boundary is greater than or equal to 1
%                             
                            [InterpCBar, Coeff_Rho] = obj.Find_CBar(xi_bar, C_Gamma, X_Cross, XSize, poly, CenterID, ka, kd, poly(CenterID).D, "X"); 
                           out4 = X_Cross;
                        elseif max(dX1,dX2)>max(dY1,dY2) %If there are more points in Y OR if there are equal points and the slope of the normal is less than or equal to 1
%                             
                            [InterpCBar, Coeff_Rho] = obj.Find_CBar(xi_bar, C_Gamma, Y_Cross, YSize, poly, CenterID, ka, kd, poly(CenterID).D, "Y");
                           out4 = Y_Cross;
                        else
                            "No case exists for nonempty X_Cross and Y_Cross"
                            InterpCBar = [];
                            Coeff_Rho = 0;
                            out4 = [0,0];
                        end
                    else 
                        "C_i, C_ii, and C_iii do not lie in the domain";
                        InterpCBar = [];
                        Coeff_Rho = 0;
                        out4 = [0,0];
                    end
            end
            out1 = InterpCBar; 
            out2 = Coeff_Rho;
            out3 = C_Gamma;
            out5 = LineID;
        end
        
        %% Create the network for each subdomain and identify which discretization points line in each channel.
        function [out] = SolveNetwork(obj, k, poly, Mtrix ,SolNetStruct, SolField, Case, Case2)

            if obj.NumChannels == 1
               switch Case
                   case "ZeroDirichlet"
                       
                       SolNet = SolNetStruct.Line;
                       
                       Size = obj.NetMtrix(1).Size - 2;

                       u0 = zeros(Size,1);
                       
                       [CGamma]=obj.FindCGamma(poly,Mtrix,SolNet, SolField, Case, Case2); %Solve for the concentration in the field at the discretization points of the line
                       
                       rho_source1 = -obj.NetMtrix(1).Sum_kd; 
                       rho_source2 = CGamma(:,1)+CGamma(:,2); 

                       check_zero = (abs(rho_source1) < 1e-10); 

                       for i=1:length(SolNet) %Solve Source term over half timestep
                         if (check_zero)
                           SolNet(i) = rho_source2(i)*k/2 + SolNet(i); 
                         else
                           SolNet(i) = (SolNet(i) + rho_source2(i)./rho_source1).* ...
                           exp(rho_source1*k/2) - (rho_source2(i)./rho_source1); 
                         end 
                       end
                       
                       bL=Mtrix.RLHS*SolNet;
                   
                       u2=Mtrix.LHS\bL; %%Solve Ax=b problem 
                       
                       for i=1:length(u2) %Solve Source term over half timestep
                          if (check_zero)
                            u0(i) = rho_source2(i)*((k/2)) + u2(i); 
                          else
                            u0(i) = (u2(i) + rho_source2(i)./rho_source1).* ...
                            exp(rho_source1*((k/2))) - (rho_source2(i)./rho_source1); 
                          end 
                        end
                       SolNetStruct.Line = u0;                     
                       out=SolNetStruct;
               end
               
            elseif obj.NumChannels == 2
                switch Case
                   case "ZeroDirichlet"
                    %%%   STOPPED
                    SolNet = SolNetStruct.Line;
                    
                    Junction = SolNetStruct.Junction;
                    
                    Size =  length(Mtrix.RhoID_X);
                    u0 = zeros(Size,1);
                    rho_source1 = zeros(Size,1);
                    
                    for i = 1 : obj.NumChannels
                        LocCGamma = find(Mtrix.ID ==i); %Find the location in the vector where the matrix ID is for the correct channel
                        rho_source1(LocCGamma) = rho_source1(LocCGamma)-obj.NetMtrix(i).Sum_kd;
                        rho_source1_Junct = obj.NetMtrix(i).Sum_kd;
                    end  
                    
                    [CGamma]=obj.FindCGamma(poly,Mtrix,SolNet, SolField, Case, Case2);
%                     [CGamma_Junct]=obj.FindCGamma(poly,Mtrix,Junct, SolField, Case, Case2);
                    rho_source2 = CGamma(:,1)+CGamma(:,2);
%                     rho_source2_Junct = CGamma_Junct(:,1)+CGamma_Junct(:,2);
                    
                   check_zero = (abs(rho_source1) < 1e-10); 
%                    check_zero_Junct = (abs(rho_source1_Junct) < 1e-10); 

                   for i=1:length(SolNet) %Solve Source term over half timestep
                     if (check_zero)
                       SolNet(i) = rho_source2(i)*k/2 + SolNet(i); 
                     else
                       SolNet(i) = (SolNet(i) + rho_source2(i)./rho_source1(i)).* ...
                       exp(rho_source1(i)*k/2) - (rho_source2(i)./rho_source1(i)); 
                     end 
                   end
                   
%                    if (check_zero) %Solve Source term over half timestep for the junction
%                        Junction = rho_source2_Junct*k/2 + Junction; 
%                    else
%                        Junction = (Junction + rho_source2_Junct./rho_source1_Junct).* ...
%                        exp(rho_source1_Junct*k/2) - (rho_source2_Junct./rho_source1_Junct); 
%                    end 

                   bL=Mtrix.RLHS*SolNet;

                   u2=Mtrix.LHS\bL; %%Solve Ax=b problem 
                   
                   

                   for i=1:length(u2) %Solve Source term over half timestep
                      if (check_zero)
                        u0(i) = rho_source2(i)*((k/2)) + u2(i); 
                      else
                        u0(i) = (u2(i) + rho_source2(i)./rho_source1(i)).* ...
                        exp(rho_source1(i)*((k/2))) - (rho_source2(i)./rho_source1(i)); 
                      end 
                   end
                    
                   Junction = Mtrix.Junction'*u0;

                   SolNetStruct.Line = u0;   
                   SolNetStruct.Junction = Junction;  
                   out=SolNetStruct;

                end       
            end
        end
        
        %% Solve for the concentration in the fields
        
        function [out1] = SolveField (obj, step, dt, poly, NMtrix, FieldMtrix, SolNetStruct, SolField, Case)
            %Given each point in C_Gamma, solve for the concentration in rho
            
            %Update Solution
            if obj.NumChannels == 1
                SolNet = [0;SolNetStruct.Line;0];
            else
                %SolNet = SolNetStruct.Line;
                Junction = SolNetStruct.Junction;
            end
            
           for i = 1:obj.NumDomains
               S = size(SolField(i).Field);
               Rho = zeros(S(1),1);
               
               for j = 1:obj.NumChannels

                   Disc = obj.NetMtrix(j).Disc;
                   DiscIn2D_X = obj.NetMtrix(j).DiscIn2D_X;
                   DiscIn2D_Y = obj.NetMtrix(j).DiscIn2D_Y;
                   
                 %  try
                        LocC_Gamma = find(FieldMtrix(i).CGamma(:,end) == j); %Locate the CGammas associated with each 
                        
                        LocRho = find(NMtrix.ID == j); %Find the points in rho that fall on this line
               
                        % ClosestPoint = knnsearch([NetMtrix.RhoID_X,NetMtrix.RhoID_X],[FieldMtrix(i).CGamma(:,2),FieldMtrix(i).CGamma(:,3)])
              
                        for k = 1:length(LocC_Gamma)
                     
                            dist = min(abs(DiscIn2D_X - FieldMtrix(i).CGamma(LocC_Gamma(k),2))); %Find the distance between the minimum point and CGamma
                    
                            Loc = find((abs(DiscIn2D_X - FieldMtrix(i).CGamma(LocC_Gamma(k),2)) - dist)<obj.tol); %Find the location of the minimum point 
                 
                            if DiscIn2D_X(Loc) - FieldMtrix(i).CGamma(LocC_Gamma(k),2)>0
                                Loc = Loc-1;
                                try
                                Gamma = Disc(Loc)+sqrt((DiscIn2D_X(Loc) - FieldMtrix(i).CGamma(LocC_Gamma(k),2))^2+(DiscIn2D_Y(Loc) - FieldMtrix(i).CGamma(LocC_Gamma(k),3))^2);
                                catch
                                    1
                                end
                            else
                                Gamma = Disc(Loc)+sqrt((DiscIn2D_X(Loc) - FieldMtrix(i).CGamma(LocC_Gamma(k),2))^2+(DiscIn2D_Y(Loc) - FieldMtrix(i).CGamma(LocC_Gamma(k),3))^2);
                            end
                            
                    if obj.NumChannels ==1
                            if Loc == 1  
                                Location = [Loc;Loc+1;Loc+2;Loc+3];       
                            elseif Loc == 2 
                                Location = [Loc-1;Loc;Loc+1;Loc+2];
                            elseif Loc == length(LocRho)-2
                                Location = [Loc-3;Loc-2;Loc-1;Loc];
                            elseif Loc == length(LocRho)-3
                                 Location = [Loc-2;Loc-1;Loc;Loc+1];
                            else
                                Location = [Loc-1;Loc;Loc+1;Loc+2]; 
                            end
                    
                            Y_Vals = Disc(Location);   %Find the values of rho at the correct locations in the discretization           

                            Coeff2=obj.LagrangeIP(Y_Vals,Gamma);%Find Lagrange Polynomials

                            rhoUpdate = sum(Coeff2'.*SolNet(Location)); %This part is messed up
                    elseif obj.NumChannels == 2
                        
                        if obj.NetMtrix(j).LocVertex == 1 %Since SolNet has no boundaries, they need to be included in the correct locations
                            SolNet = [Junction;SolNetStruct.Line(NMtrix.ID == j);0]; %Include junction first
                        else
                            SolNet = [0;SolNetStruct.Line(NMtrix.ID == j);Junction]; %Include junction last
                        end
                        
                        if Loc == 1  %Change interplolation depending on where the boundaries are
                            Location = [Loc;Loc+1;Loc+2;Loc+3];       
                        elseif Loc == 2 
                            Location = [Loc-1;Loc;Loc+1;Loc+2];
                        elseif Loc == length(SolNet)
                            Location = [Loc-3;Loc-2;Loc-1;Loc];
                        elseif Loc == length(SolNet)-1
                             Location = [Loc-2;Loc-1;Loc;Loc+1];
                        else
                            Location = [Loc-1;Loc;Loc+1;Loc+2]; 
                        end
                            

                               
                            Sol = SolNet(Location);
                            Y_Vals = Disc(Location); 
                            Coeff2=obj.LagrangeIP(Y_Vals,Gamma);%Find Lagrange Polynomials

                            rhoUpdate = sum(Coeff2'.*Sol); %This part is messed up 
                    end
                            Rho(FieldMtrix(i).CGamma(LocC_Gamma(k),1)) = Rho(FieldMtrix(i).CGamma(LocC_Gamma(k),1)) + FieldMtrix(i).CGamma(LocC_Gamma(k),4)*rhoUpdate;
                            
                            %Rho(FieldMtrix(i).CGamma(LocC_Gamma(k),1)) = Rho(FieldMtrix(i).CGamma(LocC_Gamma(k),1)) + rhoUpdate; 
                        end
                  % catch
                  %     "Insufficient facts always invite danger. Review SolveField function."
                   %    pause
                  % end
               end
               
               switch Case
                   case "LaplaceError"
                       u1 = (obj.TrueDM5(poly(i).InXNB,poly(i).InYNB) -  ((2/dt)*FieldMtrix(i).Laplace*(SolField(i).Field) + (2/dt)*Rho)); %Convergence of Laplacian
                   case "Solution"
                       if i ==1
                           [~, KPP_vals] = ode45(@(t,u) obj.KPP1(t,u), [step*dt,step*dt+(dt/2)],SolField(1).Field);
                           u2 = KPP_vals(end,:)';
                           
                           u0 = FieldMtrix(1).LHS\(FieldMtrix(1).RLHS*(u2) + 2*Rho);
                           
                           [~, KPP_vals] = ode45(@(t,u) obj.KPP1(t,u), [step*dt+(dt/2),step*dt+dt],u0);     %Solve Reaction Term over half timestep
                           u1 = KPP_vals(end,:)'; 
                       elseif i ==2
                           [~, KPP_vals] = ode45(@(t,u) obj.KPP2(t,u), [step*dt,step*dt+(dt/2)],SolField(2).Field);
                           u2 = KPP_vals(end,:)';
                           
                           u0 = FieldMtrix(2).LHS\(FieldMtrix(2).RLHS*(u2) + 2*Rho);
                           
                           [~, KPP_vals] = ode45(@(t,u) obj.KPP2(t,u), [step*dt+(dt/2),step*dt+dt],u0);     %Solve Reaction Term over half timestep
                           u1 = KPP_vals(end,:)'; 
                       end
                       
                       
               end
             
             
             
            
            %  u0 = (obj.TrueDM5(poly(i).InXNB,poly(i).InYNB));
            % u0 = 2*FieldMtrix(i).Laplace*(SolField(i).Field) + 2*Rho
             
             SolField(i).Field = u1;
           end
               
           out1 = SolField; 
            
        end
 %%       
        function [out] = FindCGamma(obj,poly,Mtrix,SolNet, SolField, Case, Case2)
            
            if obj.NumChannels == 1
               switch Case
                   case "ZeroDirichlet"
                       
                       InX1 = poly(1).InXNB;      %Set the coordinates of the points in domain 1
                       InY1 = poly(1).InYNB;
                       Sol1 = SolField(1).Field(:);
                       
                       InX2 = poly(2).InXNB;      %Set the coordinates of the points in domain 2
                       InY2 = poly(2).InYNB;
                       Sol2 = SolField(2).Field(:);
                       
                       Cross = [Mtrix.RhoID_X, Mtrix.RhoID_Y]; %Set the points of interest
                       
                       kaPointer = obj.NetMtrix.Pointer_ka; %Get the ka value for each of the surrounding subdomains
                       Loc1 = find(kaPointer(:,1)==1); %Find the location of ka for domain 1
                       Loc2 = find(kaPointer(:,1)==2); %Find the location of ka for domain 2
                       
                       Size = obj.NetMtrix(1).Size - 2;
                       
                       
                       C1a = zeros(Size,1);
                       C2a = zeros(Size,1);
                       C1 = zeros(Size,1);
                       C2 = zeros(Size,1);
                       
                       for i = 1:Size 
                           
                           [EBMCoeff1, Coeff_Rho1, C_Gamma1, C_i1, LineID1] = obj.EBM(poly, 1, Cross(i,:), Cross(i,:));
                           [EBMCoeff2, Coeff_Rho2, C_Gamma2, C_i2, LineID2] = obj.EBM(poly, 2, Cross(i,:), Cross(i,:));
                           switch Case2
                                case "Plot" %Plot results of the EBM
                                    obj.PlotEBM(poly, 1, Cross(i,:),Cross(i,:), C_Gamma1, C_i1, EBMCoeff1)
                                    obj.PlotEBM(poly, 2, Cross(i,:),Cross(i,:), C_Gamma2, C_i2, EBMCoeff2)
                           end                           

                           if ~isempty(EBMCoeff1) %Check that EBM is not empty in domain 1
                               SizeEBM = size(EBMCoeff1); %Find the number of points used to interpolate C_Bar
                              % InterpRho = [InterpRho;i,C_Gamma1(1), C_Gamma1(2), Coeff_Rho1, LineID1]; %Stores the [location in matrix, Location of C_Gamma in X and Y, Coefficient for that point]
                               for l = 1:SizeEBM(1)
                                    Loc = find(abs(EBMCoeff1(l,1)-InX1)<obj.tol &  abs(EBMCoeff1(l,2)-InY1)<obj.tol);   %Find the location of each neighbor in the vector of all points
                                    C1a(i) = C1a(i) + kaPointer(Loc1,2)*EBMCoeff1(l,3)*Sol1(Loc);                                                        %Add the correct coefficent in that spot of the matrix
                               end
                               
                           else
                               "Empty C1";
                           end
                           
                           if ~isempty(EBMCoeff2) %Check that EBM is not empty in domain 2
                               SizeEBM = size(EBMCoeff2); %Find the number of points used to interpolate C_Bar
                              % InterpRho = [InterpRho;i,C_Gamma2(1), C_Gamma2(2), Coeff_Rho2, LineID2]; %Stores the [location in matrix, Location of C_Gamma in X and Y, Coefficient for that point]
                               for l = 1:SizeEBM(1)
                                    Loc = find(abs(EBMCoeff2(l,1)-InX2)<obj.tol &  abs(EBMCoeff2(l,2)-InY2)<obj.tol);   %Find the location of each neighbor in the vector of all points
                                    C2a(i) = C2a(i) + kaPointer(Loc2,2)*EBMCoeff2(l,3)*Sol2(Loc);                                                           %Add the correct coefficent in that spot of the matrix
                               end
                               
                           else
                               "Empty C1";
                           end
                           
                           C1(i)=C1a(i)+Coeff_Rho1*SolNet(i);
                           C2(i)=C2a(i)+Coeff_Rho2*SolNet(i);
                          
                       end
               end
            elseif obj.NumChannels > 1
                switch Case
                   case "ZeroDirichlet"
                       
                       C1 = zeros(length(Mtrix.RhoID_X),1);
                       C2 = zeros(length(Mtrix.RhoID_X),1);
                       
                       for i = 1:obj.NumChannels
                           LocCGamma = find(Mtrix.ID ==i); %Find the location in the vector where the matrix ID is for the correct channel
                           Cross = [Mtrix.RhoID_X(LocCGamma), Mtrix.RhoID_Y(LocCGamma)]; %Set the points of interest
                           LocBoundaries = find(obj.Network.LineID == i); %Find the pointer to the bounding domains for this edge
                           BndDms = obj.Network.Boundary(LocBoundaries(1),:);
                           
                           InX1 = poly(BndDms(1,1)).InXNB;      %Set the coordinates of the points in domain 1
                           InY1 = poly(BndDms(1,1)).InYNB;
                           Sol1 = SolField(BndDms(1,1)).Field(:);
                       
                           InX2 = poly(BndDms(1,2)).InXNB;      %Set the coordinates of the points in domain 2
                           InY2 = poly(BndDms(1,2)).InYNB;
                           Sol2 = SolField(BndDms(1,2)).Field(:);
                           
                           SolNetTemp = SolNet(LocCGamma); %Find the values of the 
                           
                           kaPointer = obj.NetMtrix.Pointer_ka; %Get the ka value for each of the surrounding subdomains
                           Loc1 = find(kaPointer(:,1)==BndDms(1,1)); %Find the location of ka for domain 1
                           Loc2 = find(kaPointer(:,1)==BndDms(1,2)); %Find the location of ka for domain 2
                       
                           Size = length(LocCGamma);
                           
                           C1a = zeros(Size,1);
                           C2a = zeros(Size,1);
                           C1f = zeros(Size,1);
                           C2f = zeros(Size,1);
                       
                       for j = 1:Size 
                           
                           [EBMCoeff1, Coeff_Rho1, C_Gamma1, C_i1, LineID1] = obj.EBM(poly, BndDms(1,1), Cross(j,:), Cross(j,:));
                           [EBMCoeff2, Coeff_Rho2, C_Gamma2, C_i2, LineID2] = obj.EBM(poly, BndDms(1,2), Cross(j,:), Cross(j,:));
                           switch Case2
                                case "Plot" %Plot results of the EBM
                                    obj.PlotEBM(poly, 1, Cross(j,:),Cross(j,:), C_Gamma1, C_i1, EBMCoeff1)
                                    obj.PlotEBM(poly, 2, Cross(j,:),Cross(j,:), C_Gamma2, C_i2, EBMCoeff2)
                           end                           

                           if ~isempty(EBMCoeff1) %Check that EBM is not empty in domain 1
                               SizeEBM = size(EBMCoeff1); %Find the number of points used to interpolate C_Bar
                              % InterpRho = [InterpRho;i,C_Gamma1(1), C_Gamma1(2), Coeff_Rho1, LineID1]; %Stores the [location in matrix, Location of C_Gamma in X and Y, Coefficient for that point]
                               for l = 1:SizeEBM(1)
                                    Loc = find(abs(EBMCoeff1(l,1)-InX1)<obj.tol &  abs(EBMCoeff1(l,2)-InY1)<obj.tol);   %Find the location of each neighbor in the vector of all points
                                    C1a(j) = C1a(j) + kaPointer(Loc1,2)*EBMCoeff1(l,3)*Sol1(Loc);                                                        %Add the correct coefficent in that spot of the matrix
                               end
                               
                           else
                               "Empty C1";
                           end
                           
                           if ~isempty(EBMCoeff2) %Check that EBM is not empty in domain 2
                               SizeEBM = size(EBMCoeff2); %Find the number of points used to interpolate C_Bar
                              % InterpRho = [InterpRho;i,C_Gamma2(1), C_Gamma2(2), Coeff_Rho2, LineID2]; %Stores the [location in matrix, Location of C_Gamma in X and Y, Coefficient for that point]
                               for l = 1:SizeEBM(1)
                                    Loc = find(abs(EBMCoeff2(l,1)-InX2)<obj.tol &  abs(EBMCoeff2(l,2)-InY2)<obj.tol);   %Find the location of each neighbor in the vector of all points
                                    C2a(j) = C2a(j) + kaPointer(Loc2,2)*EBMCoeff2(l,3)*Sol2(Loc);                                                           %Add the correct coefficent in that spot of the matrix
                               end
                               
                           else
                               "Empty C2";
                           end
                           
                           C1f(j)=C1a(j)+Coeff_Rho1*SolNetTemp(j);
                           C2f(j)=C2a(j)+Coeff_Rho2*SolNetTemp(j);
                          
                       end
                       C1(LocCGamma) = C1(LocCGamma) +  C1f;  
                       C2(LocCGamma) = C2(LocCGamma) +  C2f;
                       end
                end
                
                
            end
            out=[C1,C2];
        end
    end
end

