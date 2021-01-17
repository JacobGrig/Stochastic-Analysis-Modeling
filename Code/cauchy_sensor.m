function res_matr = cauchy_sensor(x_0, gamma, n_rows, n_cols)

if is_not_real(x_0)
    error('x_0 should be real! Be careful!');
end

if is_not_positive(gamma)
    error('Lambda should be positive! Be careful!');
end

if is_not_natural(n_rows)
    error('The number of rows should be a natural number! Be careful!');
end

if is_not_natural(n_cols)
    error('The number of columns should be a natural number! Be careful!');
end

res_matr = x_0 + gamma * tan(pi * (rand(n_rows, n_cols) - 0.5));

end