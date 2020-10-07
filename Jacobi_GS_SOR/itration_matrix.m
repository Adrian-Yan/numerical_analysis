clear;
format long;
data;

I = eye(n);
Jacobi = I - inv(diag(diag(a)))*a;
disp(Jacobi);
GS = inv(diag(diag(a))+tril(a))*(-triu(a));
disp(GS);