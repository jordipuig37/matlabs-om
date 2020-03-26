
function [xk,dk,alk,iWk,betak,Hk] = om_uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu)
    if isd == 3
        [xk,dk,alk,iWk,Hk] = bfgs(x, f, g, epsG, kmax, c1, c2, almax, almin, rho, iW);
        betak = zeros(1, length(xk));
    elseif isd == 2
        [xk,dk,alk,iWk,betak] = conjugate_grad(x, f, g, epsG, kmax, icg, irc, c1, c2, almax, almin, rho, iW, nu);
        Hk = zeros(1, length(xk));
    elseif isd == 1
        [xk,dk,alk,iWk] = gradient(x, f, g, epsG, kmax, c1, c2, almax, almin, rho, iW);
        betak = zeros(1, length(xk));
        Hk = zeros(1, length(xk));
    else
        fprintf("isd must be in [1,3]");
    end
    function [xk,dk,alk,iWk,Hk] = bfgs(x, f, g, epsG, kmax, c1, c2, almax, almin, rho, iW)
        xk = []; dk= []; alk = [];iWk = []; Hk = [];
        I = eye(length(x));
        H = I; k = 0;
        %x_aux = x;
        while norm(g(x)) > epsG && k < kmax
            d = - H * g(x);
            dk = [dk, d];
            [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
            alk = [alk, al]; iWk = [iWk, iWout];
            x_aux = x;
            x = x_aux + al * d;
            xk = [xk, x];
            s = x - x_aux; y = g(x) - g(x_aux); p = 1/(y'*s);
            Hk = [Hk, H];
            H = (I - p*s*y') * H * (I - p*y*s') + p*s*s';
            k = k+1;
        end
    end
    function [xk,dk,alk,iWk,betak] = conjugate_grad(x, f, g, epsG, kmax, icg, irc, c1, c2, almax, almin, rho, iW, nu)
        xk = []; dk= []; alk = [];iWk = []; betak = [];
        d = - g(x); k = 0; x_ aux = x;
        while norm(g(x)) > epsG && k < kmax
            [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
            alk = [alk, al]; iWk = [iWk, iWout];
            x_aux = x;
            x = x_aux + al * d;
            xk = [xk, x];
            b = 0;
            if not(irc == 1 && mod(k, size(x, 2)) == 0) || not(irc == 2 && g(x).'* g(x_1) / norm(g(x)) >= nu)
                if icg == 1
                    b = (g(x)' * g(x)) / norm(g(x_1)); % FR
                elseif icg == 2
                    b= (g(x)' * (g(x) - g(x_1))) / norm(g(x_1));
                    b = max(0, b); % PR+
                end
            else
                b = 0;
            end
            betak = [betak, b];
            dk = [dk, d];
            d = - g(x) + b * d_1;
            k = k+1;
        end
    end
    function [xk,dk,alk,iWk] = gradient(x, f, g, epsG, kmax, c1, c2, almax, almin, rho, iW)
        xk = []; dk= []; alk = [];iWk = [];
        xk = [x];  k = 0;
        while norm(g(x)) > epsG && k < kmax
            d = -g(x);
            dk = [dk, d];
            [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
            alk = [alk, al]; iWk = [iWk, iWout];
            x = x + al*d;
            xk = [xk,x];
        end
    end
end
