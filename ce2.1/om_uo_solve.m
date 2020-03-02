function [d, b] = descent_dir(isd, x, x_1, f, g, h)
    if isd == 1
        d = - g(x);
    elseif isd == 2
        b = 0;
        if icg == 1
            b = (g(x)' * g(x)) / norm(g(x_1)); % FR
        elseif icg == 2
            b= (g(x)' * (g(x) - g(x_1))) / norm(g(x_1));
            b = max(0, b); % PR+
        end
        d = - g(x) + b * d_1; % trobem la nova direcció de descens

    elseif isd == 3
        d = 0; % encara no sé
    end
end

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

% funció principal:
function [xk,dk,alk,iWk,betak,Hk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu)
    xk = [x]; dk = []; alk = [];
    iWk = []; betak = []; Hk = [h(x)];
    k = 0;
    d_1 = 0; x_1 = 0;
    while norm(g(x)) >= 10^-6 && k < kmax
        % descent direction: %
        [d, beta] = descent_dir(isd, x, x_1, f, g, h);
        d_1 = d; % actualitzem d_1 per la proxima iteració
        % line search: %
        if iW == 0
            al = 0; % formula
        else
            [al, iWi] = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW);
        end
        % computem el pas %
        x = x + al * d;
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al];
        iWk = [iWk, iWi]; betak = [betak, beta]; Hk = [Hk, h(x)];
        k = k + 1;
    end
end
