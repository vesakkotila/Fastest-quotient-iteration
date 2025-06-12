function y = applinv(x,A,B,C,D)
%INVCH Apply the inverse using either pivoted Cholesky decomposition or
%fully permuted LU-decomposition.
%   Assuming P'*A*P=R'*R is the pivoted Cholesky decomposition of A,
%   compute y=inv(A)*x.
%   p is the permutation of vector 1:n corresponding to permutation matrix P
%   pt is the permutation of vector 1:n corresponding to permutation matrix P'

if nargin==3
    R=A;
    p=B;
    pt(p)=1:length(p);
    y=R\(R'\x(p,:));
    y=y(pt,:);
elseif nargin==5
    L=A;
    U=B;
    p=C;
    q=D;
    pt(p)=1:length(p);
    qt(q)=1:length(q);
    y=U\(L\x(p));
    y=y(qt);
else 
    error('Incorrect use of function applinv, see description.')
end