function dJ=dleg(n,t)
% DLEGANDRE: find dJ={L'_j(t_i)}=[L'_0 L'_1 ... L'_n](t). 
% size(dJ)=[numel(t),n+1].
% L_k = Legendre polinomial of order k on [-1,1],
% L'_n=0.5(n+1)P^(1,1)_{n-1} -- Jcobi polinomial.

t=t(:);
one =ones(size(t));
if n==1, dJ=one; return; end

dJ=zeros(numel(t),n);  % find dJ_k=P^(1,1)_k
dJ(:,1)=one; dJ(:,2)=2*t;
for k=1:n-2
  dJ(:,k+2)=(((2*k+3)*(k+2)/(k+1))*t.*dJ(:,k+1)-(k+2)*dJ(:,k))/(k+3);
end
N=(1:n)'; d=0.5*(N+1);
dJ=[zeros(numel(t),1) dJ*diag(d)];
