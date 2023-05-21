function [result] = unstatic_right(t, Ux, A, R, L)
n = length(Ux)/4;
result = A*Ux-[applied_force(t, n, R, L)];
% tnc = (R*L)^(-1);
%  = [tnc*temp(1:n);tnc*temp(n+1:2*n);tnc*temp(2*n+1:3*n);tnc*temp(3*n+1:4*n)];
end