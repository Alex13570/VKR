function J=leg(n,t)
% LEGENDRE: find J={L_j(t_i)}=[L_0 L_1 ... L_n](t), 
% size(J)=[numel(t),n+1].
% L_k = Legendre polinomial of order k on [-1,1].

t=t(:);
one=ones(size(t));
if n==0, J=one; return; end

J=zeros(numel(t),n+1);
J(:,1)=one; J(:,2)=t;
for k=1:n-1
  J(:,k+2)=((2*k+1)*t.*J(:,k+1)-k*J(:,k))/(k+1);
end
