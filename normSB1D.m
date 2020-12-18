function [out] = normSB1D(A, dx, Case)
% Grid norms as defined in Leveque A.5 pg 252 
switch Case
    case 2
        out=(dx*sum(abs(A).^2))^(1/2);
    case 0
        out=max(abs(A));
end
end