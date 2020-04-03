function [xk, dk, alk, iWk, Hk, tauk] = modified_newton(isd, x, f, g, h, epsG, kmax, almax,almin,rho,c1,c2,iW, delta)
    xk = [x]; dk= []; alk = [];iWk = []; Hk = []; tauk = [];
    k = 0;
    while norm(g(x)) > epsG && k < kmax
        B = 0;
        if isd == 5
            B = SD(x, h, delta);
            Hk = cat(3, Hk, B);
        else
            [B, tau] = cholo(x, h);
            Hk = cat(3, Hk, B);
            tauk = [tauk, tau];
        end
        d = - inv(B) * g(x);
        dk = [dk,d];
        [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
        alk = [alk, al]; iWk = [iWk, iWout];
        x = x + al * d;
        xk = [xk, x];
        k = k+1;
    end
    if isd == 5
        tauk = zeros(1, length(xk));
    end

    function B = SD(x, h, delta)
        % Q' A Q = h(x)
        [Q, A] = schur(h(x));
        n = length(x);
        for i = 1:n
            Abar(i,i) = max(delta, A(i,i));
        end
        B = Q' * Abar * Q;
    end
    function [B, tau] = cholo(x, h)
        upBound = norm(h(x), 'fro');
        tau = 0; k = 0;
        while not_exists_rr(h(x), tau)
            k = k+1;
            tau = (1.01 - 1/(2^k)) * upBound;
        end
        B = h(x) + tau * eye(length(x));
    end
    function boolea = not_exists_rr(h, tau)
        [r, boolea] = chol(h + eye(length(h(:,1))) * tau);
    end
end
