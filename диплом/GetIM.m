function [mIM1,mIM2,mIM3,x]=GetIM(cN,vX)
%%
n=cN;
a=vX(1);
b=vX(end);
[t,w] = qgauss(n,-1,1);
x=(a+b)/2+(b-a)/2*t;
J =leg(n-1,t);
DJ=dleg(n-1,t);
[L,R, d] = img(a,b,t,w,J(:,2:end),DJ(:,2:end));
mIM1 = L;
mIM2 = R;
global vALF0
vALF0 = x;
mIM3 = d;
%--------------------------------------------------------------------------
