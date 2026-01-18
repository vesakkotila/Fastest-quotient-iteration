% Figures for Example 8
figure
ll=eigs(M,N,3,0);
ln=length(quotient1{1});
semilogy(0:(ln-1),quotient1{1}-ll(1),'b--')
hold
semilogy((ln-1):(ln-2+length(quotient2{1})),abs([quotient1{1}(end) quotient2{1}(2:end)]-ll(1)),'b')

title('''Error'' in eigenvalue approximation (compared to eigs)')
xlabel('Iteration')
legend('\lambda_1: PDM','\lambda_1: OQI')

figure
semilogy(0:(ln-1),quality1{1},'b--')
hold
semilogy((ln-1):(ln-2+length(quality2{1})),[quality1{1}(end) quality2{1}(2:end)],'b')
title('Quality of the eigenvector in terms of \sigma_2([w_1 w_2])')
xlabel('Iteration')
legend('\lambda_1: PDM','\lambda_1: OQI')