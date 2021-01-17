function t = QS_mod_mod(lambda_0, T)

t = 0;

t_cur = 0;

while true
    xi = exprnd(1 / (2 * lambda_0));
    t_cur = t_cur + xi;
    if t_cur > T
        break
    end
    eta = bern_gen((1 + cos(t_cur)) / 2, 1, 1);
    
    if (eta)
        t = [t, t_cur];
    end
end