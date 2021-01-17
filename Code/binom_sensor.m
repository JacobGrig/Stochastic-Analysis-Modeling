function res_vect = binom_sensor(p, n_par, quant)

if is_not_prob(p)
    error('Probability should be a real number between zero and one! Be careful!');
end

if is_not_natural(n_par)
    error('The second parameter of binomial distribution should be a natural number! Be careful!');
end

if is_not_natural(quant)
    error('The quantity of selection should be a natural number! Be careful!');
end

bern_matr = bern_gen(p, n_par, quant);
res_vect = sum(bern_matr, 1)';
end