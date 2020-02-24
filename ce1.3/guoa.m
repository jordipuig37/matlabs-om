% [start] guoa algorithm
% donats un x inicial, una funcio a minimitzar f, i els paràmetres de guoa
% retorna la x que minimitza trobada. parametre max_it per el màxim d'iteracions
% i parametre newton per si es vol seguir el metode del gradient (0) o newton (1)
function [x, k] = guoa(x, f, g, gg, almax, almin, c1, c2, rho, max_it, newton)
    k = 0; iW = 1;
    while norm(g(x)) >= 10^-6 && k < max_it
        % gradient method:
        d = - g(x);
        if newton == 1
        % newton method
            d = d / gg(x);
        end
        al = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW);
        x=x+al*d; k=k+1;
    end
end
