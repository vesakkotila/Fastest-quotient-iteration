% Figures for Example 9
figure
ll=eigs(M,N,3,0)
ln=length(quotient1{1});
semilogy(0:(ln-1),quotient1{1}-ll(1),'b--')
hold
semilogy((ln-1):(ln-2+length(quotient2{1})),abs([quotient1{1}(end) quotient2{1}(2:end)]-ll(1)),'b')
ln2=length(quotient1{2});
semilogy(0:(ln2-1),abs(quotient1{2}-ll(2)),'r--')
semilogy((ln2-1):(ln2-2+length(quotient2{2})),abs([quotient1{2}(end) quotient2{2}(2:end)]-ll(2)),'r')

title('''Error'' in eigenvalue approximation (compared to eigs)')
xlabel('Iteration')
legend('\lambda_1: PDM','\lambda_1: OQI','\lambda_2: PDM','\lambda_2: OQI')

figure
semilogy(0:(ln-1),quality1{1},'b--')
hold
semilogy((ln-1):(ln-2+length(quality2{1})),[quality1{1}(end) quality2{1}(2:end)],'b')
semilogy(0:(ln2-1),quality1{2},'r--')
semilogy((ln2-1):(ln2-2+length(quality2{2})),[quality1{2}(end) quality2{2}(2:end)],'r')
title('Quality of the eigenvector in terms of \sigma_2([w_1 w_2])')
xlabel('Iteration')
legend('\lambda_1: PDM','\lambda_1: OQI','\lambda_2: PDM','\lambda_2: OQI')