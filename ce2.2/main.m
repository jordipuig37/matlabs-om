clear;
% Problem
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

% Input parameters.
 % Stopping criterium:
epsG = 10^-6; kmax = 1500;
 % Linesearch:
almax = 2; almin = 10^-3; rho=0.5; c1=0.01; c2=0.45; iW = 2;
 % Search direction:
isd = 1; icg = 2; irc = 0 ; nu = 0.1;

% Optimization
[xk,dk,alk,iWk,betak,Hk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu);
save('uo_FDM_CE21.mat','f','g','h','epsG','kmax','almax','almin','rho','c1','c2','iW','isd','icg','irc','nu','xk','dk','alk','iWk','betak');

% Output
niter = size(xk,2); xo = xk(:,niter); fk = []; gk = []; rk = []; gdk = [];
for k = 1:niter
    fk = [fk,f(xk(:,k))]; gk=[gk,g(xk(:,k))]; % f(xk) and g(xk)
end
for k = 1:niter-1
    rk = [rk,(fk(k+1)-f(xo))/(fk(k)-f(xo))]; % Rate of convergence
end
for k = 1:niter-1
    gdk = [gdk,gk(:,k)'*dk(:,k)]; % Descent condition
end
rk=[rk,NaN]; gdk = [gdk,0];
fprintf('[om_uo_FDM_CE21]\n');
fprintf('epsG= %3.1e, kmax= %4d\n', epsG,kmax);
fprintf('almax= %2d, almin= %3.1e, rho= %4.2f\n',almax,almin,rho);
fprintf('c1= %3.2f, c2= %3.2f, iW= %1d\n',c1,c2,iW);
fprintf('isd= %1d, icg= %1d, irc= %1d, nu= %3.1f\n',isd,icg,irc,nu);
fprintf('    k   x(1)     x(2)   iW   g''*d      ||g||   r\n');
rk(niter) = 0;
for k = 1:niter-1
    fprintf('%5d %7.4f %7.4f %3d %+3.1e %4.2e %3.1e %3.1e\n', k, xk(1,k), xk(2,k), iWk(k), gdk(k), norm(gk(:,k)), rk(k), rk(k)'*rk(k+1));
end
fprintf('   k x(1)  x(2)    iW  g''*d   ||g||   r\n[om_uo_FDM_CE21]\n');
xylim=[0 0 0 0];
subplot(1,2,1);
om_uo_solve_plot(f, xk, gk, xylim, 1, 0);
subplot(1,2,2);
om_uo_solve_plot(f, xk, gk, xylim, 2, 0);
