%This script computes the first two eigenvalues for the waveguide problem of
%Example 9. It is assumed that the following files are in the MATLAB path:
%  applinv.m, mgs.m and oqiP.m along with Zmodel.m, which produces the matrices.
% First, four preconditioned descent steps are taken to produce the 
% quotients in vector quotient1, after which Algorithm 3 (for the smallest 
% eigenvalue) or Algorithm 1 (for other eigenvalues) is invoked to have
% iterates in quotient2.

Zmodel;
n=size(M,1);

clear qual* quot* x

K=4; %number of descent steps
Numeig=2; %number of eigenvalues to be computed
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
%innerP=@(v,w) w'*applP(v);


% 1st eigenvector
q=rand(n,1);
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
        x=x/norm(x);
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
x=x/norm(x);

%2nd,... eigenvector
for l=2:Numeig
    % Set mu equal to the previous eigenvalue, add a small number to
    % prevent stagnation.
    mu=quotient2{l-1}(end)+1e-7;
    
    [R,FLAG,p]=chol(M-mu*N,'vector'); % for applying Z=inv(M-mu*N)
    
    if FLAG==0 % if positive definite,     
        applZ= @(q) applinv(q,R,p); % use Cholesky
    else
        [L,U,p,r]=lu(M-mu*N,'vector'); % use fully permuted LU
        applZ=@(q) applinv(q,L,U,p,r);
    end
    
    q=rand(n,1);
    q=q-x(:,1:(l-1))*(x(:,1:(l-1))'*q);
    q=q/normP(q);
    x(:,l)=applZ(q);
    x(:,l)=x(:,l)/norm(x(:,l));
    quotient1{l}(1)= normP(M*x(:,l))/normP(N*x(:,l));
    w1=M*x(:,l);
    w1=w1/normP(w1);
    w2=N*x(:,l);
    w2=w2/normP(w2);
    quality1{l}(1)=svds([w1 w2],1,"smallest");

    

   
    Mhatq=M*applZ(q); %Mhat*q=M*Z*q, Z=inv(M-mu*N)
    Nhatq=N*applZ(q); 
    for k=1:K
        %Compute the conjugate co-gradient
        %a=N\q-2*mu*(M\q)+mu^2*M\(N*(M\q));
        a=applZ(Mhatq); %Mhat'*P*Mhatq=Z*Mhatq
        b=applZ(N*applP(Nhatq));   %Nhat'*P*Nhatq=Z*N*P*Nhatq=inv(M-muN)*N*P*Nhatq
        beta=(q'*a)/(q'*b);
        g=a-(q'*a)/(q'*b)*b;
        % Orthogonalize against the previous eigenvectors
        g=g-x(:,1:(l-1))*(x(:,1:(l-1))'*g);
        
        %Solve preconditioned eigenvalueproblem in the basis {q,qhat}
        % with (q,qhat)=0. Orthogonalize qhat against the previous
        % eigenvectors
        if normP(g)>1e-15
            qhat=mgs(g,q);
            qhat=qhat-x(:,1:(l-1))*(x(:,1:(l-1))'*qhat);
            qhat=qhat/normP(qhat);
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
            q=q-x(:,1:(l-1))*(x(:,1:(l-1))'*q);
            q=q/normP(q);
            Mhatq=M*applZ(q);
            Nhatq=N*applZ(q);
            x(:,l)=applZ(q);
            x(:,l)=x(:,l)/norm(x(:,l));
            % The approximate quotient for the lth eigenvalue,
            % after k preconditioned descent steps.
            quotient1{l}(k+1)= normP(M*x(:,l))/normP(N*x(:,l));
            % Determine the quality of the approximation by 
            % studying the linear dependece of Mx and Nx.
            w1=M*x(:,l);
            w1=w1/normP(w1);
            w2=N*x(:,l);
            w2=w2/normP(w2);
            quality1{l}(k+1)=svds([w1 w2],1,"smallest");
        else 
            break
        end 
        
       
    end 

    % Switch to Algorithm 1
    [quotient2{l},x(:,l),quality2{l}]=oqiP(M,N,27000,2,x(:,l),1e-10);
    x(:,l)=x(:,l)/norm(x(:,l));
end 
    
   