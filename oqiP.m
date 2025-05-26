function [l,q,quality]=oqiP(M,N,mu,ipm,q0,tol)

% OQIP Optimal quotient iteration with respect to norm P.
% 
% Approximate the eigenvalue of (M-mu*N)q=l*N by 
% Algorithm 1 in M. Huhtanen, V. Kotila, P. Uusitalo: "FASTEST 
% QUOTIENT ITERATION WITH VARIATIONAL PRINCIPLES FOR SELF-ADJOINT 
% EIGENVALUE PROBLEMS" .
% 
% INPUT: 
% ipm       inner product matrix: choose 1 for P=inv(M), 2 for P=inv(N)
% q0        starting vector (optional)
% tol       tolerance for linear independence of w1=Mq/||Mq|| and 
%           w2=Nq/||Nq||, measured by svd([w1 w2])
% OUTPUT:
% l     approximate eigenvalue
% q     corresponding approximate eigenvector
% quality   quality of the eigenvector measured by svd([w1 w2]),
% 


switch nargin
    case 5
        tol=1e-10;
    case 4
        q0=rand(size(M,2),1);
        tol=1e-10;
    case 3
        ipm=1;
        q0=rand(size(M,2),1);
        tol=1e-10;        
end

q=q0;
M=M-mu*N;

Mq=M*q;
Nq=N*q;
 if ipm==1
    [R,FLAG,p]=chol(M,'vector'); % for applying Z=inv(M) and normP
    if FLAG==0
        applP=@(q) applinv(q,R,p);
    else
        error('N not positive definite!')
    end
   
    nMq=sqrt(q'*Mq); %||Mq||_P
    nNq=normP(Nq,R,p); %||Nq||_P
    w1=Mq/nMq;
    w2=Nq/nNq;
    w1w2=(q'*Nq)/(nMq*nNq);
    S=svd([Mq/norm(Mq) Nq/norm(Nq)]); %S(2) is the measure for lin. dependence of Mq and Nq
    k=0;
    l=w1w2/abs(w1w2)*nMq/nNq;
    while S(2)>tol
        k=k+1;
        z=1/sqrt(2+2*abs(w1w2))*(w1w2/abs(w1w2)*w1+w2);
        qhat=(M-l(k)*N)\z;
        q=qhat/normP(qhat,R,p);
        if k>100
            break;
        end
        Mq=M*q;
        Nq=N*q;
        nMq=sqrt(q'*Mq);
        nNq=normP(Nq,R,p);
        w1=Mq/nMq;
        w2=Nq/nNq;
        w1w2=(q'*Nq)/(nMq*nNq);
        S=svd([Mq/norm(Mq) Nq/norm(Nq)]);
        quality(k)=S(2);
        l(k+1)=w1w2/abs(w1w2)*nMq/nNq;
    end

 end

if ipm==2
    [R,FLAG,p]=chol(N,'vector'); % for applying Z=inv(M) and normP
    if FLAG==0
        applP=@(q) applinv(q,R,p);
    else
        error('N not positive definite!')
    end
    %normP=@(q) sqrt(q'*applP(q));
    nNq=normP(Nq,R,p);
    nMq=normP(Mq,R,p);
    w1=Mq/nMq;
    w2=Nq/nNq;
    w1w2=(q'*Mq)/(nMq*nNq);
    l=w1w2/abs(w1w2)*nMq/nNq;
    sigma2=svds([Mq/norm(Mq) Nq/norm(Nq)],1,"smallest");
    sigmaold=2;
    k=0;
    while (sigma2>tol) && (sigma2<sigmaold)
        k=k+1;
        z=1/sqrt(2+2*abs(w1w2))*(w1w2/abs(w1w2)*w1+w2);
        qhat=(M-l(k)*N)\z;
        q=qhat/normP(qhat,R,p);
        if k>100
            break;
        end
        Mq=M*q;
        Nq=N*q;
        nMq=normP(Mq,R,p);
        nNq=normP(Nq,R,p);
        w1=Mq/nMq;
        w2=Nq/nNq;
        w1w2=(q'*Mq)/(nMq*nNq);
        sigmaold=sigma2;
        sigma2=svds([Mq/norm(Mq) Nq/norm(Nq)],1,"smallest");
        quality(k)=sigma2;
        l(k+1)=w1w2/abs(w1w2)*nMq/nNq;
    end
    l=l+mu;
end



 end