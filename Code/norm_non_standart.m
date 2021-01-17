function res_vect = norm_non_standart(mu, sigma2, quant)

if is_not_positive(sigma2)
    error("Sigma^2 should be positive! Be careful!");
end

if is_not_real(mu)
    error("Mu should be real! Be careful!");
end

if is_not_natural(quant)
    error("An amount of elements should be natural! Be careful");
end

if ~rem(quant, 2)
    [xi_vect, eta_vect] = norm_pair(quant / 2);
    res_vect = [xi_vect, eta_vect];
else
    [xi_vect, eta_vect] = norm_pair((quant + 1) / 2);
    res_vect = [xi_vect, eta_vect(1 : (end - 1))];
end
res_vect = sqrt(sigma2) * res_vect + mu;
end