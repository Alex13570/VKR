function [d,c] = quadg(n,a,b,dgts)
% QGAUSS: find nodes and weights of n-point 
% Gauss quadrature rule on [a,b].  
% d = nodes=zeroes of Legendre polinomial L_n(t) shifted on [a,b], 
% c = weights. 
% d,c  are column vectors

digits(dgts)

if nargin==1, a=-1; b=1; end
if n==1, d=(a+b)/2; c=2; return; end

j=vpa(1:n-1);
beta=vpa(j./sqrt(4*j.^2-1));

A=vpa(diag(beta,-1)+diag(beta,1));

% Find nodes x and weights w for [-1,1].
[V,D]=eig(vpa(A));
d=vpa(a+0.5*(b-a)*(diag(D)+1));

[d,k]=sort(double(d)); % nodes on [a,b]

w=double(vpa(2*V(1,k).^2));   % weights on [-1,1]
c=double(vpa(0.5*(b-a)*w'));  % weights on [a,b]

