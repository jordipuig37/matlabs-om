n=2;
rng(123456);
Q=rand(n);
[V,D] = eig(triu(Q)+triu(Q)');
Q=V*diag(diag(max(D,1)))*V';
b=rand(n,1); xo = Q^-1*b; 
x = rand(n,1);
f = @(x) x'*Q*x/2-b'*x;
g = @(x) Q*x-b;
h = @(x) Q;
