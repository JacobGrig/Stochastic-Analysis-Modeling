function res_vect = pois_binom_sensor(lambda, N, quant)

p = lambda / N;

if is_not_prob(p)
    error('Probability seems to be not a real number between zero and one! Be careful!');
end

res_vect = binom_sensor(p, N, quant);

end