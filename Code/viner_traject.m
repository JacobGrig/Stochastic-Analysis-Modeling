function [t_vect, W_vect] = viner_traject(alpha, n)

t_vect = nan(1, 2 ^ n + 1);
W_vect = nan(1, 2 ^ n + 1);
used_indices = false(1, 2 ^ n + 1);

t_vect(1) = 0;
W_vect(1) = 0;

t_vect(end) = 1;
W_vect(end) = randn(1, 1);

used_indices(1) = true;
used_indices(end) = true;

counter = 1;

for i = 1 : n
    
    t_cur = t_vect(used_indices);
    W_cur = W_vect(used_indices);
    
    t_new = alpha .* t_cur(2 : end) + (1 - alpha) .* t_cur(1 : (end - 1));
    
    means = alpha .* W_cur(2 : end) + (1 - alpha) .* W_cur(1 : (end - 1));
    vars = alpha .* (1 - alpha) .* (t_cur(2 : end) - t_cur(1 : (end - 1)));
    
    W_new = sqrt(vars) .* randn(1, counter) + means;
    
    counter = counter * 2;
   
    old_indices = find(used_indices);
    new_indices = (old_indices(1 : (end - 1)) + old_indices(2 : end)) / 2;
    
    used_indices(new_indices) = true;
    
    t_vect(new_indices) = t_new;
    W_vect(new_indices) = W_new;
end