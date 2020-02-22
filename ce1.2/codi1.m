clear;
Q = [4 0; 0 1];
f = @(x) (1/2)*x'*Q*x; g = @(x) Q*x;
x = [1;2]; k = 1; iW =2;
almax = 1; almin = 0.01; rho = 0.5;
c1 = 0.1; c2 = 0.5;

while norm(g(x)) >= 10^-6 && k <= 100
    d = -g(x);
    al = -((Q*x)' * d)/(d'*Q*d);
    % al = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
    fprintf('%3d %10.6f %10.6f %10.6f \n', k, x(1), x(2),al);
    x=x+al*d; k=k+1;
end

fprintf('%3d %10.6f %10.6f\n', k, x(1), x(2));

%%% no hi ha millora, al contrari, tardem 24 operacions i tot i aixi
%%% no arribem al (0,0). Aixi demostra per contra exemple que la
%%% alpha optima no sempre dona solucions i execucions �ptimes.
