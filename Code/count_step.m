function res_vect = count_step(probs)

len = size(probs, 2);

prob_sums = cumsum(probs);

x_vect = rand(1, len);

[res_vect, ~] = find([x_vect < prob_sums(1, :); ...
    (x_vect < prob_sums(2 : end, :)) & (x_vect >= prob_sums(1 : (end - 1), :))]);

end