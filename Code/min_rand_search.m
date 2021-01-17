function [min_elem, min_value, funcs] = min_rand_search(func, number)
radiuses = sqrt(rand(1, number));
alphas = 2 * pi * rand(1, number);

samples = [radiuses .* cos(alphas); radiuses .* sin(alphas)];
funcs = func(samples);

[min_value, min_ind] = min(funcs);
min_elem = samples(:, min_ind);
end