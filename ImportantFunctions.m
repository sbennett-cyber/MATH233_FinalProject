classdef ImportantFunctions
   properties

   end
   methods
       function obj = ImportantFunctions()
             
       end
       %%
       function [out] = normSB1D(obj, A, dx, Case)
            % Grid norms as defined in Leveque A.5 pg 252 
            switch Case
                case 2
                    out=(dx*sum(abs(A).^2))^(1/2);
                case 0
                    out=max(abs(A));
            end
       end
       %%
       function [out] = normSB(obj,A, dx, dy, Case)
        % Grid norms as defined in Leveque A.5 pg 252 
        switch Case
            case 2
                out=(dx*dy*sum(sum(abs(A).^2)))^(1/2);
            case 0
                out=max(max(abs(A)));
        end
        end
       %%
       function [out, out1, out2] = ErrorBtwnSol(obj, R1, R2, R3, R4, R5, dx, dt, Case1)
        %%% Shayna Bennett - 12/13/2020
        %%% This function finds the order of a method by comparing the solutions at
        %%% the final time when either the spatial step or the timestep are varied. 
        %%% This follows page 257 of LeVeque's "Finite Difference Methods for
        %%% Ordinary and Partial Differential Equations".
        Len = length(dx);
            switch Case1
                case "Space"
                    RR2=R1-R2(1:2:end);
                    RR3=R2-R3(1:2:end);
                    RR4=R3-R4(1:2:end);
                    RR5=R4-R5(1:2:end);
                case "Time"
                    RR2=R1-R2;
                    RR3=R2-R3;
                    RR4=R3-R4;
                    RR5=R4-R5;
            end
                Error0=[0,obj.normSB1D(RR2,dx(1),0),obj.normSB1D(RR3,dx(2),0),obj.normSB1D(RR4,dx(3),0),obj.normSB1D(RR5,dx(4),0)]';
                Error2=[0,obj.normSB1D(RR2,dx(1),2),obj.normSB1D(RR3,dx(2),2),obj.normSB1D(RR4,dx(3),2),obj.normSB1D(RR5,dx(4),2)]';
                Order0=zeros(Len,1);
                Order2=zeros(Len,1);
                for a=3:(Len)
                    Order0(a,1)=log2((Error0(a-1))/(Error0(a)));
                    Order2(a,1)=log2((Error2(a-1))/(Error2(a)));
                end
                out = table(dx,dt, Error0, Order0, Error2, Order2);
                out1=Order0;
                out2=Order2;       
        end
%%
   function [out, out1, out2] = ErrorBtwnSol2D(obj, R1, R2, R3, R4, R5, dx, dt, Case1)
        %%% Shayna Bennett - 12/13/2020
        %%% This function finds the order of a method by comparing the solutions at
        %%% the final time when either the spatial step or the timestep are varied. 
        %%% This follows page 257 of LeVeque's "Finite Difference Methods for
        %%% Ordinary and Partial Differential Equations".
        Len = length(dx);
            switch Case1
                case "Space"
                    RR2=R1-R2(1:2:end,1:2:end);
                    RR3=R2-R3(1:2:end,1:2:end);
                    RR4=R3-R4(1:2:end,1:2:end);
                    RR5=R4-R5(1:2:end,1:2:end);
                case "Time"
                    RR2=R1-R2;
                    RR3=R2-R3;
                    RR4=R3-R4;
                    RR5=R4-R5;
            end
                Error0=[0,obj.normSB(RR2,dx(1),dx(1),0),obj.normSB(RR3,dx(2),dx(2),0),obj.normSB(RR4,dx(3),dx(3),0),obj.normSB(RR5,dx(4),dx(4),0)]';
                Error2=[0,obj.normSB(RR2,dx(1),dx(1),2),obj.normSB(RR3,dx(2),dx(2),2),obj.normSB(RR4,dx(3),dx(3),2),obj.normSB(RR5,dx(4),dx(4),2)]';
                Order0=zeros(Len,1);
                Order2=zeros(Len,1);
                for a=3:(Len)
                    Order0(a,1)=log2((Error0(a-1))/(Error0(a)));
                    Order2(a,1)=log2((Error2(a-1))/(Error2(a)));
                end
                out = table(dx,dt, Error0, Order0, Error2, Order2);
                out1=Order0;
                out2=Order2;       
        end
%%
   function [out, out1, out2] = ErrorWTrueSol2D(obj, R1, R2, R3, R4, R5, dx, dt)
        %%% Shayna Bennett - 12/13/2020
        %%% This function finds the order of a method by comparing the solutions at
        %%% the final time when either the spatial step or the timestep are varied. 
        %%% This follows page 257 of LeVeque's "Finite Difference Methods for
        %%% Ordinary and Partial Differential Equations".
        Len = length(dx);
                Error0=[obj.normSB(R1,dx(1),dx(1),0),obj.normSB(R2,dx(2),dx(2),0),obj.normSB(R3,dx(3),dx(3),0),obj.normSB(R4,dx(4),dx(4),0),obj.normSB(R5,dx(5),dx(5),0)]';
                Error2=[obj.normSB(R1,dx(1),dx(1),2),obj.normSB(R2,dx(2),dx(2),2),obj.normSB(R3,dx(3),dx(3),2),obj.normSB(R4,dx(4),dx(4),2),obj.normSB(R5,dx(5),dx(5),2)]';
                Order0=zeros(Len,1);
                Order2=zeros(Len,1);
                for a=2:(Len)
                    Order0(a,1)=log2((Error0(a-1))/(Error0(a)));
                    Order2(a,1)=log2((Error2(a-1))/(Error2(a)));
                end
                out = table(dx,dt, Error0, Order0, Error2, Order2);
                out1=Order0;
                out2=Order2;       
        end

function [out, out1, out2] = ErrorWTrueSol1D(obj, R1, R2, R3, R4, R5, dx, dt)
        %%% Shayna Bennett - 12/13/2020
        %%% This function finds the order of a method by comparing the solutions at
        %%% the final time when either the spatial step or the timestep are varied. 
        %%% This follows page 257 of LeVeque's "Finite Difference Methods for
        %%% Ordinary and Partial Differential Equations".
        Len = length(dx);
                Error0=[obj.normSB1D(R1,dx(1),0),obj.normSB1D(R2,dx(2),0),obj.normSB1D(R3,dx(3),0),obj.normSB1D(R4,dx(4),0),obj.normSB1D(R5,dx(5),0)]';
                Error2=[obj.normSB1D(R1,dx(1),2),obj.normSB1D(R2,dx(2),2),obj.normSB1D(R3,dx(3),2),obj.normSB1D(R4,dx(4),2),obj.normSB1D(R5,dx(5),2)]';
                Order0=zeros(Len,1);
                Order2=zeros(Len,1);
                for a=2:(Len)
                    Order0(a,1)=log2((Error0(a-1))/(Error0(a)));
                    Order2(a,1)=log2((Error2(a-1))/(Error2(a)));
                end
                out = table(dx,dt, Error0, Order0, Error2, Order2);
                out1=Order0;
                out2=Order2;       
        end
   end
end