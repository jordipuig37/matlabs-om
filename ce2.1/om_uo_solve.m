% funció principal:
function [xk,dk,alk,iWk,betak,Hk]= om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu)
    xk = [x]; dk = []; alk = [];
    k = 0;
    d_1 = 0;
    while norm(g(x)) >= 10^-6 && k < kmax
        % descent direction: %
        if isd == 1
            d = - g(x);

        elseif isd == 2
            if icg == 1
                b = % FR
            elseif icg == 2
                b = % PR+
            end
            d = - g(x) + b * d_1; % trobem la nova direcció de descens
            d_1 = d; % actualitzem d_1 per la proxima iteració

        elseif isd == 3
            d = 0; % encara no sé
        end
        % line search: %
        if iW == 0
            al = ; % formula
        else
            al = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW);
        end
        x = x + al * d;
        xk = [xk, x]; dk = [dk, d]; alk = [alk, al]
        k = k + 1;
    end
end
