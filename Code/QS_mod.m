function t = QS_mod(lambda_0, T)

% fsolve_opts = optimoptions('fsolve', 'Display','none');

Lambda = @(t, y) lambda_0 .* (t + sin(t)) - y;
Lambda_t = @(t) Lambda(t, 0);
Lambda_m1 = @(y, t0) fzero(@(t) Lambda(t, y), t0);

t = 0;

while true
    E = exprnd(1, 1, 1);
    t_new = Lambda_m1(E + Lambda_t(t(end)), t(end));
    if (t_new < t(end))
        disp('lol');
    end
    if t_new > T
        break
    end
    t = [t, t_new];
end

end