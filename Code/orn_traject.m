function [t_vect, W_vect] = orn_traject(sigma, lambda, n, T)

t_vect = nan(1, 2 ^ n + 1);
W_vect = nan(1, 2 ^ n + 1);
used_indices = false(1, 2 ^ n + 1);

t_vect(1) = 0;
W_vect(1) = sigma * randn(1, 1);

t_vect(end) = T;
W_vect(end) = W_vect(1) * exp(-lambda * T) + ...
    sigma * sqrt(1 - exp(-2 * lambda * T)) * randn(1, 1);

used_indices(1) = true;
used_indices(end) = true;

counter = 1;

for i = 1 : n
    
    t_cur = t_vect(used_indices);
    W_cur = W_vect(used_indices);
    
    x1s = W_cur(1 : (end - 1));
    x2s = W_cur(2 : end);
    
    t1s = t_cur(1 : (end - 1));
    t2s = t_cur(2 : end);
    
    t_new = (t2s + t1s) ./ 2;
    
    means = (x1s + x2s) .* (exp(-lambda .* (t2s - t1s) ./ 2)) ./ ...
        (1 + exp(-lambda .* (t2s - t1s)));
    vars = sigma .^ 2 .* (1 - exp(-lambda .* (t2s - t1s))) ./ ...
        (1 + exp(-lambda .* (t2s - t1s)));
    
    W_new = sqrt(vars) .* randn(1, counter) + means;
    
    counter = counter * 2;
   
    old_indices = find(used_indices);
    new_indices = (old_indices(1 : (end - 1)) + old_indices(2 : end)) / 2;
    
    used_indices(new_indices) = true;
    
    t_vect(new_indices) = t_new;
    W_vect(new_indices) = W_new;
end
end