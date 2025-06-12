function [q,h] = mgs(a,Q)
%MGS Modified Gram-Schmidt
%   Given a vector a and orthonormal column vectors Q(:,1),...,Q(:,n) of Q,
%   calculate the scalar projections h(1),...,h(n) and the new orthogonal
%   vector q=v/h(n+1), where v=a-h(1)Q(:,1)-...-h(n)Q(:,n) and
%   h(n+1)=norm(v).
%   
[m,n]=size(Q);
v=a;
for j=1:n,
    h(j,1)=dot(v,Q(:,j));
    v=v-h(j,1)*Q(:,j);
end
h(n+1,1)=norm(v);

if abs(h(n+1,1))>1e-15,
    q=v/h(n+1,1);
else
    q=v;
    warning('Nearly linearly dependent vectors')
end
