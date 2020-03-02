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
