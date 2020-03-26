% funció principal:
function [xk,dk,alk,iWk,betak,Hk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu)
    xk = [x]; dk = []; alk = [];
    iWk = []; betak = []; Hk = [h(x)];
    k = 0; H =eye();
    d_1 = 0; x_1 = 0;
    while norm(g(x)) >= epsG && k < kmax
        % descent direction: %
        [d, b, H] = descent_dir(isd, irc, x, x_1, f, g, H, k);
        d_1 = d; % actualitzem d_1 per la proxima iteració
        % line search: %
        if iW == 0
            al = 0; % formula
            iWi = 3;
        else
            [al, iWi] = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW);
        end
        % computem el pas %
        x_1 = x;
        x = x + al * d;
        % posem tota la informació al dia %
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al];
        iWk = [iWk, iWi]; betak = [betak, b]; Hk = [Hk, h(x)];
        k = k + 1;
    end
end
