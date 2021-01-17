function res_vect = match_pen(n)

if is_not_natural(n)
    error('The second parameter of binomial distribution should be a natural number! Be careful!');
end

bern_vect = 2 * bern_gen(0.5, n, 1) - 1;
res_vect = cumsum(bern_vect) / sqrt(n);
end