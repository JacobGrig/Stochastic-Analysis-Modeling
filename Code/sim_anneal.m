function [min_elem, min_value, probs, counter] = sim_anneal(func, m, n, D, T_inits, x0, limits, T_limit, epsilon)

c = m .* exp(-n ./ D);

T_func = @(k) T_inits .* exp(-c .* k .^ D);

i = 1;
x_cur = x0;
Ti = T_func(1);

probs = [];

% delta = 1e-1;

counter = 0;

while true
    counter = counter + 1;
    g_old = func(x_cur);
    alphas = rand(1, D).';
    z = sign(alphas - 0.5) .* Ti(2 : end) .* ((1 + 1 ./ Ti(2 : end)) .^ abs(2 .* alphas - 1) - 1);
    
    if sum(isinf(z))
        break;
    end
    
    x_cand = x_cur + z .* (limits(:, 2) - limits(:, 1));
    
    if sum((x_cand < limits(:, 1)) | (x_cand > limits(:, 2)))
%         if min(abs(x_cur - x_cand)) < delta
        left_border = (x_cand < limits(:, 1));
        right_border = (x_cand > limits(:, 2));
        x_cand(left_border) = limits(left_border, 1);
        x_cand(right_border) = limits(right_border, 2);
%         else
%             continue;
%         end
    end
    
    g_new = func(x_cand);
    
    differ = g_new - g_old;
    
    if (differ < 0)
        
        probs = [probs, 1];
        if (norm(x_cur - x_cand) < epsilon) && (sum(Ti < T_limit))
            x_cur = x_cand;
%             disp('break by x-diff');
            break;
        end
        x_cur = x_cand;
        i = i + 1;
        Ti = T_func(i);
        if sum(Ti < T_limit)
            break
        end
    else
        p = exp(-differ / Ti(1));
        probs = [probs, p];
        if logical(bern_gen(p, 1, 1))
%             if (norm(x_cur - x_cand) < epsilon) && (sum(Ti < T_limit))
%                 break;
%             end
            x_cur = x_cand;
            i = i + 1;
            Ti = T_func(i);
%             if sum(Ti < T_limit)
%                 break
%             end
        end
    end
end
min_elem = x_cur;
min_value = g_new;
end