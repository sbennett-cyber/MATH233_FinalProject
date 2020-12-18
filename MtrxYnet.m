function [out] = MtrxYnet(alpha,NNx,Channel)
%This function returns the centered difference discretization in X and can
%be used for either a line or a plane.
%%Takes as input a coefficent alpha, spatial step, and number of points in x and y directions

%%% Assume Zero Dirichlet BC on all outer channels %%%

A = sparse(diag(ones(NNx,1)*(-2*alpha))) + ...
    sparse(diag(ones(NNx-1,1)*alpha,1)) + ...
    sparse(diag(ones(NNx-1,1)*alpha,-1));
%A(end,end-1:end)=[8/9, -14/9]               %Corrects for interfacial point
%A2(1,1:2)=[8/9, -14/9] 
u_dm_diag = []; 
%u_Side=[];
%u_side=[zeros(NNx-1,1);1];
%u_bottom=[zeros(NNx-2,1);0.5;-2];
% u_B1=alpha*[zeros(NNx-2,1);8/9;-14/9;-1/9;4/9;zeros(NNx-2,1);4/9;-1/9;zeros(NNx-2,1)]';
% u_B2=alpha*[zeros(NNx-2,1);-1/9;4/9;zeros(NNx-2,1);8/9;-14/9;4/9;-1/9;zeros(NNx-2,1)]';
% u_B3=alpha*[zeros(NNx-2,1);-1/9;4/9;-1/9;4/9;zeros(NNx-2,1);-14/9;8/9;zeros(NNx-2,1)]';

u_fix=alpha*(2/9)*[zeros(NNx-2,1);-1/2;2;zeros(NNx-2,1);-1/2;2;2;-1/2;zeros(NNx-2,1)]';

for i=1:Channel
  u_dm_diag = blkdiag(u_dm_diag,A);
end
u_dm_diag(NNx,:) = u_dm_diag(NNx,:)+u_fix;
u_dm_diag(2*NNx,:) = u_dm_diag(2*NNx,:)+u_fix;
u_dm_diag(2*NNx+1,:) = u_dm_diag(2*NNx+1,:)+u_fix;

out=u_dm_diag;
end
