function res_matr = exp_sensor(lambda, n_rows, n_cols)

if is_not_positive(lambda)
    error('Lambda should be positive! Be careful!');
end

if is_not_natural(n_rows)
    error('The number of rows should be a natural number! Be careful!');
end

if is_not_natural(n_cols)
    error('The number of columns should be a natural number! Be careful!');
end

res_matr = -log(rand(n_rows, n_cols)) ./ lambda;

end