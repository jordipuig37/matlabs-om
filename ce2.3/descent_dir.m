function [d, b, H] = descent_dir(isd, irc, x, x_1, d_1, f, g, H, nu, k)
    b = 0;
    if isd == 1
        % gradient method
        d = - g(x);
    elseif isd == 2
        % conjugate gradient (if not restart-> get beta, else beta=0)
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
        d = - g(x) + b * d_1; % trobem la nova direcció de descens

    elseif isd == 3
        % quasi-newton
        y = g(x) - g(x_1);
        s = x - x_1;
        p = 1 / (y.' * s);
        H  = (I - p*s*y.') * H * (I - p*s*y.') +  p*s*y.';
        d = - g(x)/H;
    else
        d = 0;
    end
end
