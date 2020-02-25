%%% Classe 20/02/2020 %%%
%%% implementacio CE1.1: %%%
% definim la funcio i els paràmetres:
Q = [4 0; 0 1]; f = @(x) (1/2)*x'*Q*x; g = @(x) Q*x;
x = [1;2]; d = [-4;-2];
almax= 1.0; almin= 10^-6; rho = 0.5; c1 = 0.1; c2 = 0.5; iW=1;

[al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);



% [start] Alg. BLS 
function [al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)
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
    if swc(al)
        iWout = 3;
    elseif wc(al)
        iWout = 2;
    elseif wc1(al)
        iWout = 1;
    end
end
% [end] Alg. BLS 

