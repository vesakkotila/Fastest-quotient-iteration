function [l,q,quality]=alg3P(M,N,ipm,q0,tol)

% Approximate the smallest eigenvalue with Algorithm 3.
%tol=1e-10;
% alg3P Optimal quotient iteration for eigenvector approximation associated
%with the smallest eigenvalue of a positive definite self-adjoint 
% eigenvalue problem.  
% 
% Algorithm 3 in M. Huhtanen, V. Kotila, P. Uusitalo: "FASTEST 
% QUOTIENT ITERATION WITH VARIATIONAL PRINCIPLES FOR SELF-ADJOINT 
% EIGENVALUE PROBLEMS" .
% 
% INPUT: 
% ipm       inner product matrix: choose 1 for P=inv(M), 2 for P=inv(N)
% q0        starting vector (optional)
% tol       tolerance for linear independence of w1=Mq/||Mq|| and 
%           w2=Nq/||Nq||, measured by svd([w1 w2]) (optional)
% OUTPUT:
% l         iteration history for the approximate eigenvalue 
% q         corresponding approximate eigenvector (only the last)
% quality   iteration history for the quality of the eigenvector,
%           measured by svd([w1 w2]),
% 

q=q0;
Mq=M*q;
Nq=N*q;
 if ipm==1 %P=inv(M)
    [R,FLAG,p]=chol(M,'vector'); % for applying normP
    if FLAG==0
        normP=@(q) sqrt(q'*applinv(q,R,p));
    else
        error('M not positive definite!')
    end
    nMq=sqrt(q'*Mq);
    nNq=normP(Nq);
    w1=Mq/nMq;
    w2=Nq/nNq;
    w1w2=(q'*Nq)/(nMq*nNq);
    S=svd([Mq/norm(Mq) Nq/norm(Nq)]);
    k=0;
    l=1/Algo2(N,M,R,p,q);   % Initial value for the quotient. 
    quality1{1}(1)=S(2);    % Initial value for the quality of the 
                            % eigenvector
    while S(2)>tol
        k=k+1;
        z=1/sqrt(2+2*abs(w1w2))*(w1w2/abs(w1w2)*w1+w2);
        qhat=(M-l(k)*N)\z;
        q=qhat/normP(qhat);
        if k>100
            break;
        end
        Mq=M*q;
        Nq=N*q;
        nMq=sqrt(q'*Mq);
        nNq=normP(Nq);
        w1=Mq/nMq;
        w2=Nq/nNq;
        w1w2=(q'*Nq)/(nMq*nNq);
        S=svd([Mq/norm(Mq) Nq/norm(Nq)]);
        quality(k+1)=S(2);
        l(k+1)=1/Algo2(N,M,R,p,q);
    end

 end

if ipm==2 % P=inv(N)
    [R,FLAG,p]=chol(N,'vector'); % for applying Z=inv(M) and normP
    if FLAG==0
        normP=@(q) sqrt(q'*applinv(q,R,p));
    else
        error('N not positive definite!')
    end
    nNq=sqrt(q'*Nq);
    nMq=normP(Mq);
    w1=Mq/nMq;
    w2=Nq/nNq;
    w1w2=(q'*Mq)/(nMq*nNq);
    S=svd([Mq/norm(Mq) Nq/norm(Nq)]);
    k=0;
    l=1/Algo2(N,M,R,p,q);   % Initial value for the quotient. 
    quality1{1}(1)=S(2);    % Initial value for the quality of the 
                            % eigenvector                      
    while S(2)>tol
        k=k+1
        l=w1w2/abs(w1w2)*nMq/nNq;
        z=1/sqrt(2+2*abs(w1w2))*(w1w2/abs(w1w2)*w1+w2);
        qhat=(M-l*N)\z;
        q=qhat/normP(qhat);
        if k>100
            break;
        end
        Mq=M*q;
        Nq=N*q;
        nMq=sqrt(q'*Mq);
        nNq=normP(Nq);
        w1=Mq/nMq;
        w2=Nq/nNq;
        w1w2=(q'*Mq)/(nMq*nNq);
        S=svd([w1 w2])
        quality(k+1)=S(2);
        l(k+1)=1/Algo2(N,M,R,p,q);
    end

end
end