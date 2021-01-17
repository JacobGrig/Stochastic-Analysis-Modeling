function [xi_vect, eta_vect] = norm_pair(quant)

if is_not_natural(quant)
    error('The second parameter of binomial distribution should be a natural number! Be careful!');
end

unif_vect = 2 * pi * rand(1, quant);
exp_vect = exp_sensor(0.5, 1, quant);

xi_vect = sqrt(exp_vect) .* cos(unif_vect);
eta_vect = sqrt(exp_vect) .* sin(unif_vect);

end