function alpha=Algo2(M,N,R,p,q)

% Algorithm 2: find the best quotient approximating the largest eigenvalue
% of Mx=lambda*Nx given an approximate eigenvector q. Use the given Cholesky
% factorization P^(-1)=R*R' for computing the norm ||x||_P=sqrt((x,Px)).

normP=@(q) sqrt(q'*applinv(q,R,p));

v=M*q;
w=N*q;
l=normP(w);
mu=normP(v)/(2*l); % midpoint of the spectrum
muold=0;
while abs(mu-muold)>eps
    alpha=normP(v-mu*w)/l+mu;
    muold=mu;
    mu=alpha/2;
end