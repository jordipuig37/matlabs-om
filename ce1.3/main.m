% exercici 3 %
f = @(x) x^2 + 2 * sin(x);
g = @(x) 2 * cos(x) + 2 * x;
gg = @(x) 2 - 2 * sin(x);

almax = 1; almin = 0.01; rho = 0.75;
c1 = 0.1; c2 = 0.5;
% initial point
x = 10;
% 0 = gradient method, 1 = newton method
newt = 1;

[optima, it] = guoa(x, f, g, gg, almax, almin, c1, c2, rho, 10000, newt);

fprintf("Mínim %d trobat en %d operacions des de %d\n", [optima, it, x])
