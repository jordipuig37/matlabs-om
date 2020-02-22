% exercici CE1.1
function k = w1(x, d, f, g, al, c1)
    k = f(x + al * d) <= f(x) + c1 * g(x).' * d * al;
end

function k = w2(x, d, g, al, c2)
    k = g(x + al*d).' * d >= c2 * g(x).' * d;
end

function k = swc(x, d, g, al, c2)
    k = abs(g(x + al*d).' * d) >= abs(c2 * g(x).' * d);
end

function b = satisfied(iw, al, x, d, f, g, c1, c2)
    if iw == 1
        b = w1(x, d, f, g, al, c1, c2) & w2(x, d, f, g, al, c1, c2);
    elseif iw == 2
       b = w1(x, d, f, g, al, c1, c2) & swc(x, d, f, g, al, c1, c2);
    end
end

function [al,iWout] = BLS(x, d, f, g, almax, almin, rho, c1, c2, iW)
    al = almax;
    while not(satisfied(iW, al, x, d, f, g, c1, c2)) | al < almin
        al = al * rho;
        iwout = iw;
    end
end

        
Q = [4 0; 0 1];
f = @(x) (1/2)*x'*Q*x;
g = @(x) Q*x;
x = [1;2];
d = [-4;-2];
almax = 1.0;
almin = 10^-6;
rho = 0.5;
c1 = 0.1;
c2 = 0.5;
iw = 1;

BLS(x, d, f, g, almax, almin, rho, c2, c2, iw)
