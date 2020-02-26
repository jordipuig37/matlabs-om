clear;
xo=[1;1];
Q=[10, 1; 1, 10]; b=Q*xo;
Q=[10, 2; 2, 1];  b=Q*xo;
f = @(x) x'*Q*x/2-b'*x; g = @(x) Q*x-b; h = @(x) Q;
x  = [0;6]; fo= f(xo);
xk = [x]; dk= []; alk= []; fk=[f(x)]; gk=[g(x)]; rk=[]; k=1;
epsG = 10^-6;

while norm(g(x)) >= epsG
    d = -g(x); al = -(Q*x-b)'*d/(d'*Q*d); dk = [dk,d]; alk= [alk,al];
    x = x + al*d; k=k+1;
    xk = [xk,x]; fk = [fk,f(x)]; gk = [gk,g(x)]; rk = [rk,(fk(k)-fo)/(fk(k-1)-fo)];
end
rk=[rk,NaN];
fprintf('[om_uo_FDM_Ex21a_GM]\n  k    x(1)    x(2)   ||g||       r\n');
for i=1:k
    fprintf('%3d %7.4f %7.4f %3.1e  %6.4f\n', i, xk(1,i), xk(2,i), norm(gk(:,i)), rk(i));
end
la = eig(h(x)); rUB= ((la(2)-la(1))/(la(2)+la(1)))^2;
fprintf('       la=[%3.1f, %3.1f]-> rUB= %5.4f\n',la(1),la(2),rUB);
