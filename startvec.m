clear qual* quot* x

K=3; %number of descent steps
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
innerP=@(v,w) w'*applP(v);


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
    %a=N\q-2*mu*(M\q)+mu^2*M\(N*(M\q));
    a=applZ(Mhatq); %Mhat'*P*Mhatq=Z*Mhatq
    b=applZ(N*applP(Nhatq));   %Nhat'*P*Nhatq=Z*N*P*Nhatq=inv(M-muN)*N*P*Nhatq
    g=a-(q'*a)/(q'*b)*b;
    %g=g-(x'*g)*x/norm(x)^2;
    
    if normP(g)>1e-15
        qhat=mgs(g,q); %Modified Gram-Schmidt
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
        qold=q; % for testing
        Mhatq=M*applZ(q);
        Nhatq=N*applZ(q);
        x=applZ(q);
        x=x/norm(x);
    else
        break
    end

    quotient1{1}(k+1)= normP(M*x)/normP(N*x);
    w1=M*x;
    w1=w1/normP(w1);
    w2=N*x;
    w2=w2/normP(w2);
    quality1{1}(k+1)=svds([w1 w2],1,"smallest");
    
end


[quotient2{1},x,quality2{1}]=oqiP(M,N,0,1,x,1e-10);
x=x/norm(x);

for l=2:Numeig
%2nd,... eigenvector
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
        a=applZ(Mhatq); %Mhat'*P*Mhatq=Z*Mhatq
        b=applZ(N*applP(Nhatq));   %Nhat'*P*Nhatq=Z*N*P*Nhatq=inv(M-muN)*N*P*Nhatq
        beta=(q'*a)/(q'*b);
        g=a-(q'*a)/(q'*b)*b;
        g=g-x(:,1:(l-1))*(x(:,1:(l-1))'*g);
        
 
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
            quotient1{l}(k+1)= normP(M*x(:,l))/normP(N*x(:,l));
            w1=M*x(:,l);
            w1=w1/normP(w1);
            w2=N*x(:,l);
            w2=w2/normP(w2);
            quality1{l}(k+1)=svds([w1 w2],1,"smallest");
        else % else normPg
            break
        end % end if normPg
        
       
    end % end for k=

       
    [quotient2{l},x(:,l),quality2{l}]=oqiP(M,N,0,1,x(:,l),1e-10);
    x(:,l)=x(:,l)/norm(x(:,l));
end %end l=2:Numeig
    
   