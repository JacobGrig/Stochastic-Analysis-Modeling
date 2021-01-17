function [t, W] = insur(x_m, k, W_0, c, lambda, N, T)

t = 0;

while true
    exps = exprnd(1 / lambda, 1, N);
    new_t = t(end) + cumsum(exps);
    last_ind = find(new_t >= T, 1);
    if ~isempty(last_ind)
        t = [t, repelem(new_t(1 : (last_ind - 1)), 2)];
        break;
    else
        t = [t, repelem(new_t, 2)];
    end
end

number_of_falls = (length(t) - 1) / 2;

falls = x_m .* (rand(1, number_of_falls) .^ (-1 / k));

sumfalls = [0, cumsum(falls)];

W = W_0 + c .* t;

W(3 : 2 : end) = W(3 : 2 : end) - sumfalls(2 : end);
W(2 : 2 : end) = W(2 : 2 : end) - sumfalls(1 : (end - 1));

last_ind = find(W < 0, 1);

if ~isempty(last_ind)
    W = W(1 : last_ind);
    t = t(1 : last_ind);
end