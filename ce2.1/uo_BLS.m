% [start] Alg. BLS
function [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)
    % iWout = 0: al does not satisfy any WC
    % iWout = 1: al satisfies (WC1)
    % iWout = 2: al satisfies WC
    % iWout = 3: al satisfies SWC
    wc1 = @(al) f(x + al * d) <= f(x) + c1 * g(x).' * d * al;
    wc2 = @(al) g(x + al*d).' * d >= c2 * g(x).' * d;
    swc2 = @(al) abs(g(x + al*d).' * d) <= abs(c2 * g(x).' * d);

    wc = @(al) wc1(al) & wc2(al);
    swc = @(al) wc1(al) & swc2(al);

    if iW == 1
        cond = @(al) wc(al);
    elseif iW == 2
        cond = @(al) swc(al);
    end

    al = almax;
    while al > almin && not(cond(al))
        al = al * rho;
    end

    iWout = 0;
    if wc1(al)
        iWout = iWout +1;
    elseif wc2(al)
        iWout = iWout +1;
    elseif swc(al)
        iWout = iWout +1;
    end
end
% [end] Alg. BLS
