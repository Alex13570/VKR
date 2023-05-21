function [d,c,v] = qgauss(n,a,b)
% QGAUSS: find nodes d, weights v of n-point 
% Gauss quadrature rule on [a,b] and barycentric weights v.  
% d,c,v  are column vectors

if nargin==1, a=-1; b=1; end
if n==1, d=(a+b)/2; c=2; return; end

j=1:n-1;
beta=j./sqrt(4*j.^2-1);

A=diag(beta,-1)+diag(beta,1);

% Find nodes x and weights w for [-1,1].
[V,D]=eig(A);
[x,k]=sort(diag(D));
w=(2*V(1,k).^2)';

v = sqrt((1-x.^2).*w);        
v = v./max(v);
v(2:2:n) = -v(2:2:n);

c=0.5*(b-a)*w;       % --> transform to the interval [a,b]
d=a+0.5*(b-a)*(x+1);
