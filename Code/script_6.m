% Task 6

%% initials

clc;
clear;

n = 10000000;

div = 100;

N = 7;

eta = 0.95;

%% calculation
tic;

f = @(x) pi .^ 5 * (exp(-1 ./ (2 .^ 7 .* prod(x .^ 2)))) ./ (prod(x .^ 2));
%g = @(x) 1 ./ (pi .^ 5) .* exp(-sum(x .^ 2));

%f_vect = zeros(1, n);

% f_vect = f(x_vect);

summa = 0;
sum2 = 0;

results_cur = zeros(1, n / div);

for i = div : div : n
    x_vect = sqrt(0.5) .* randn(10, div);
    curs = f(x_vect);
    summa = summa + sum(curs);
    sum2 = sum2 + sum(curs .^ 2);
    if (~rem(i, 10000000))
        disp(i);
    end
    results_cur(i / div) = summa / i;
end

% results_cur = cumsum(f_vect) ./ (1 : n);

res = summa / n;

sigma2 = sum2 / n - res ^ 2;

err = norminv((eta + 1) / 2) .* sqrt(sigma2) ./ sqrt(n);

disp(['result of Monte-Carlo is ', num2str(res), ', error equals to ', num2str(err)]);

toc;
%%
% f_new = @(t) exp(-(sum(tan(pi ./ 2 .* t) .^ 2, 2) + 1 ./ (2 .^ 7 .* prod(tan(pi ./ 2 .* t) .^ 2, 2)))) ./ ...
%              (prod(tan(pi ./ 2 .* t) .^ 2, 2) .* prod(cos(pi ./ 2 .* t) .^ 2, 2));

tic;

f_new = @(t) exp(-(sum(tan(pi ./ 2 .* t) .^ 2, 2) + 1 ./ (2 .^ 7 .* prod(tan(pi ./ 2 .* t) .^ 2, 2)))) ./ ...
              (prod(sin(pi ./ 2 .* t) .^ 2, 2));
         
splitting_N = linspace(-1, 1, N + 1);
splitting_N = splitting_N(2 : end);

% splitting_t = cartesian(splitting_N, splitting_N, splitting_N, splitting_N, splitting_N, ...
%                         splitting_N, splitting_N, splitting_N, splitting_N, splitting_N) - (1 ./ (2 .* N));

splitting_t = cartesian(splitting_N, splitting_N, splitting_N, splitting_N, splitting_N) - (1 ./ (2 .* N));

res_sum = 0;

for i = 1 : length(splitting_t)
    splitting_t_mod = [splitting_t, repmat(splitting_t(i, :), N^5, 1)];
%     splitting_t_mod = [splitting_t(i, :), splitting_t(j, :)];
    f_vect_new = f_new(splitting_t_mod);
    res_sum = res_sum + sum(f_vect_new);
end
                    
factor = (pi / N) ^ 10;

% f_vect_new = f_new(splitting_t);

res_rect = factor * res_sum;

toc;

%% visualization

figure;

loglog(div : div : n, results_cur, 'b');
hold on;

loglog(div : div : n, ones(1, length(div : div : n)) * 124.83, 'r');

grid on;
xlabel('n');
ylabel('I_n');

legend('integrals', 'exact value = 124.83');

xlim([500, n]);

hold off;