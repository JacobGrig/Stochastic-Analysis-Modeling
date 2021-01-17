function [min_elem, min_value] = grad_des(func, diffs, x0, lambda_0, epsilon)

x_old = x0;
while (func(x0 - lambda_0 .* diffs(x0)) > func(x0))
    lambda_0 = lambda_0 / 2;
    if lambda_0 * norm(diffs(x0)) < epsilon
        min_elem = x0;
        min_value = func(x0);
        return;
    end
end
x_cur = x0 - lambda_0 .* diffs(x0);
lambda = lambda_0;

while (norm(x_old - x_cur) >= epsilon)
    lambda = lambda / 2;
    x_new = x_cur - lambda * diffs(x_cur);
    if (func(x_new) <= func(x_cur))
        lambda = lambda_0;
        x_old = x_cur;
        x_cur = x_new;
        continue;
    end
    if lambda * norm(diffs(x_cur)) < epsilon
        break;
    end
end
min_elem = x_cur;
min_value = func(x_cur);
end