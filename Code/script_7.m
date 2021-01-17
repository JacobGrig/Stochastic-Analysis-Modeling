% Task 7

%% initials

clc;

N = 10000;
p = 0.99;

init_T0 = 10;
init_T1 = 10;
init_T2 = 10;

T_inits = [init_T0; init_T1; init_T2];

x1_int = [-10, 10];
x2_int = [-10, 10];

limits = [x1_int; x2_int];

m0 = 1;
n0 = 1;

m1 = 1;
n1 = 1;

m2 = 1;
n2 = 1;

m = [m0; m1; m2];
n = [n0; n1; n2];

D = 2;
d = 1;

x0 = [0; 0];
T_limit = 1e-30;

f_func = @(x) func_f(x);

g_func = @(x) (x(1, :) - 1) .^ 2 + 100 .* (x(2, :) - x(1, :) .^ 2).^ 2;

dgdx1_func = @(x) 2 .* (x(1, :) - 1) + 200 .* (x(2, :) - x(1, :) .^ 2) .* (-2 .* x(1, :));
dgdx2_func = @(x) 200 .* (x(2, :) - x(1, :) .^ 2);

lambda_0 = 10;
epsilon1 = 1e-4;
epsilon2 = 1e-2;

%% calculation

[min_elem_f, min_value_f, values] = min_rand_search(f_func, N);
err = sqrt((10 + sqrt(10))^2 + 17 * 10^2) * sqrt(p / N);

%%

sim_results = zeros(1, d);
sim_res_elem = zeros(2, d);

tic;

for i = 1 : d
    try
        [min_elem_g, min_value_g, probs, counter] = sim_anneal(g_func, m, n, D, T_inits, x0, limits, T_limit, epsilon1);
    catch
        sim_results(i) = sim_results(i - 1);
    end
    sim_results(i) = min_value_g;
    sim_res_elem(:, i) = min_elem_g;
    if (~rem(i, 100))
        disp(i);
    end
end

toc;

%%
tic;
[min_elem_g_mod, min_value_g_mod] = grad_des(g_func, @(x) [dgdx1_func(x); dgdx2_func(x)], x0, lambda_0, epsilon2);
toc;
%% visualization

factor = 1;

splitting = factor : factor : N;

figure;

mins = zeros(1, N / factor);

for i = splitting
    mins(i / factor) = min(values(1 : i));
end

stairs(splitting, mins);
set(gca, 'XScale', 'log');
grid on;

xlabel('n');
ylabel('min');