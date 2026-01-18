% This script computes the smallest eigenvalue for the BCSST13-problem of
% Example 8. It is assumed that the following files are in the MATLAB path:
% applinv.m, mgs.m and algo3P.m along with mmread.m, bcsstk13.mtx and 
% bcsstm13.mtx from Matrix Market. 
% First, three preconditioned descent steps are taken to produce the 
% quotients in vector quotient1, after which Algorithm 3 is invoked to have
% iterates in quotient2.

M=mmread('bcsstk13.mtx');
N=mmread('bcsstm13.mtx');
n=size(M,1);


q0=rand(n,1);
q=q0;
clear qual* quot* x

K=3; %number of descent steps

% P=inv(M)

mu=0;

[R,FLAG,p]=chol(M-mu*N,'vector'); % for applying Z=inv(M-mu*N)
if FLAG==0 % if positive definite,     
 applZ= @(q) applinv(q,R,p); % use Cholesky
else
 [L,U,p,r]=lu(M-mu*N,'vector'); % use fully permuted LU
 applZ=@(q) applinv(q,L,U,p,r);
end
[R2,FLAG2,p2]=chol(M,'vector'); % for applying P and normP
applP=@(q) applinv(q,R2,p2);
normP=@(q) sqrt(q'*applP(q));



% 1st eigenvector

q=q/normP(q);
x=applZ(q);
x=x/norm(x);
quotient1{1}(1)= normP(M*x)/normP(N*x);
w1=M*x;
w1=w1/norm(w1);
w2=N*x;
w2=w2/normP(w2);
quality1{1}(1)=svds([w1 w2],1,"smallest");

Mhatq=M*applZ(q); %Mhat*q=M*Z*q, Z=inv(M-mu*N)
Nhatq=N*applZ(q); 
for k=1:K
    %Compute the conjugate co-gradient
    %a=N\q-2*mu*(M\q)+mu^2*M\(N*(M\q));
    a=applZ(Mhatq); %Mhat'*P*Mhatq=Z*Mhatq
    b=applZ(N*applP(Nhatq));   %Nhat'*P*Nhatq=Z*N*P*Nhatq=inv(M-muN)*N*P*Nhatq
    g=a-(q'*a)/(q'*b)*b;
    
    %Solve preconditioned eigenvalueproblem in the basis {q,qhat}
    % with (q,qhat)=0.
    if normP(g)>1e-15  
        qhat=mgs(g,q); %Use the modified Gram-Schmidt -method
        Mhatqhat=M*applZ(qhat);
        Nhatqhat=N*(applZ(qhat));
        A=[Mhatq Mhatqhat]'*applP([Mhatq Mhatqhat]);
        B=[Nhatq Nhatqhat]'*applP([Nhatq Nhatqhat]);
        [v,lam]=eig(A,B,'vector');
        if all(imag(lam)<1e-10) && all(real(lam)>0)
            [~,ind]=min(real(lam));
        else
            warning(['Av=lambda*Bv not self adjoint for k=',num2str(k),'. Using last good value.'])
            break
        end
        v=v(:,ind);
        q=[q qhat]*v;
        q=q/normP(q);
        Mhatq=M*applZ(q);
        Nhatq=N*applZ(q);
        x=applZ(q);
        x=x/normP(x);
    else
        break
    end
    
    % The approximate quotient for the 1st eigenvalue,
    % after k preconditioned descent steps.
    quotient1{1}(k+1)= normP(M*x)/normP(N*x);
    % Determine the quality of the approximation by 
    % studying the linear dependece of Mx and Nx.
    w1=M*x;
    w1=w1/normP(w1);
    w2=N*x;
    w2=w2/normP(w2);
    quality1{1}(k+1)=svds([w1 w2],1,"smallest");
    
end

% Switch to Algorithm 3
[quotient2{1},x,quality2{1}]=alg3P(M,N,1,x,1e-10);


   