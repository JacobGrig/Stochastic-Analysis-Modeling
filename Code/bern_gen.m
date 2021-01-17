function res_matr = bern_gen(p, n_rows, n_cols)

if prod(is_not_prob(p))
    error('Probability should be a real number between zero and one! Be careful!');
end

if is_not_natural(n_rows)
    error('The number of rows should be a natural number! Be careful!');
end

if is_not_natural(n_cols)
    error('The number of columns should be a natural number! Be careful!');
end

res_matr = double(rand(n_rows, n_cols) < p);
end