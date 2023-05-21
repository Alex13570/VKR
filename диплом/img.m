function [L,R, d]=img(a,b,t,d,P,dP)
% IMG: find left (L) and right (R) integrating matix on [a,b].
% L={l_ij}=int_a^{x_i} l_j(x)dx,
% R={l_ij}=int_{x_i}^b l_j(x)dx,
% where l_j(x)= Lagrange bazis function with nodes x=(a+b)/2+t*(b-a)/2,
% t,d=nodes and weights of n-point Gauss quadrature on [-1,1]
%  
% P =[L_1 ... L_{n-1}](t)-- Legendre polinomials on t
% size(P)=[numel(t),n-1]
% dP =[L'_1 ... L'_{n-1}](t), size(DP)=[numel(t),n-1]

n=numel(t);
k=1:n-1;
K=(2*k+1)./(k.*(k+1));
DK=diag(K);
W=diag(1-t.^2);
d=(0.5*(b-a))*d;

L=(1+repmat(t,1,n)-W*dP*DK*P')*diag(d/2);
R=repmat(d',n,1)-L;
%K=(repmat(t,1,n)-W*dP*DK*P')/2

