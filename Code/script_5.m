% Task 5

%% initials

clc;

mu = 1;
sigma2 = 1;

n = 100000;
N = 10000;

alpha = 0.03;

a = 1;
b = 2;

%% calculation

norm_vect = norm_non_standart(mu, sigma2, n);

splitting_n = 1 : n;
emp_norm = zeros(1, N);

for i = 1 : N
    norm_vect = norm_non_standart(mu, sigma2, n);
    emp_norm(i) = sqrt(n / sigma2) * (sum(norm_vect) / n - mu);
end
% emp_norm = sqrt(splitting_n / sigma2) .* (cumsum(norm_vect) ./ splitting_n - mu);


%%
norm_vect = norm_non_standart(mu, sigma2, n);

splitting_n = 1 : n;

means = cumsum(norm_vect) ./ (splitting_n);

f_emp_norm = (1 : N) ./ N;

vars = cumsum(norm_vect .^ 2) ./ splitting_n - means .^ 2;
s2 = (splitting_n ./ (splitting_n - 1)) .* vars;

err_1 = tinv(1 - alpha / 2, splitting_n - 1) .* sqrt(s2) ./ sqrt(splitting_n);

left_mean = means - err_1;
right_mean = means + err_1;

left_var = (splitting_n - 1) .* s2 ./ (chi2inv((1 + (1 - alpha)) / 2, splitting_n - 1));
right_var = (splitting_n - 1) .* s2 ./ (chi2inv((1 - (1 - alpha)) / 2, splitting_n - 1));

negs = s2 - left_var;
pos = right_var - s2;

disp('OK');

%%
splitting_n_mod = 1 : n;

cauchy_vect = cauchy_sensor(a, b, 1, n);

means_cauchy = cumsum(cauchy_vect) ./ (splitting_n_mod);

f_means_cauchy = (1 : N) ./ N;

cauchy_means = zeros(1, N);

for i = 1 : N
    cauchy_vect = cauchy_sensor(a, b, 1, n);
    cauchy_means(i) = sum(cauchy_vect) ./ n;
end

F_cauchy_theor = @(x) 1 / pi * atan((x - a) ./ b) + 0.5;

%% visualization

figure;

semilogx(splitting_n, means, 'r');

hold on;

semilogx(splitting_n, ones(1, n) * mu, 'b--', 'LineWidth', 1.5);

xlabel('n');
ylabel('S_n');

legend('Empirical S_n', '\mu');

axis([10, 10^5, 0.6, 1.6]);

grid on;

hold off;

figure;

plot(sort(emp_norm), f_emp_norm, 'b');

hold on;
plot(sort(emp_norm), normcdf(sort(emp_norm)), 'r--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

hold off;

figure;

plot(splitting_n(1000 : end), left_mean(1000 : end), 'b');
hold on;

plot_1 = errorbar(splitting_n(1000 : 1000 : end), means(1000 : 1000 : end), err_1(1000 : 1000 : end));

plot(splitting_n(1000 : end), right_mean(1000 : end), 'b');
plot_2 = plot(splitting_n(1000 : end), ones(1, n - 1000 + 1) * mu, 'r', 'LineWidth', 1.5);

xlabel('n');
ylabel('confidence intervals');

legend([plot_1 plot_2], {'probable means', 'exact value'});

grid on;

hold off;

figure;

plot(splitting_n(1000 : end), left_var(1000 : end), 'b');
hold on;

plot_1 = errorbar(splitting_n(1000 : 1000 : end), s2(1000 : 1000 : end), negs(1000 : 1000 : end), pos(1000 : 1000 : end));

plot(splitting_n(1000 : end), right_var(1000 : end), 'b');
plot_2 = plot(splitting_n(1000 : end), ones(1, n - 1000 + 1) * sigma2, 'r', 'LineWidth', 1.5);

xlabel('n');
ylabel('confidence intervals');

legend([plot_1 plot_2], {'probable vars', 'exact value'});

grid on;

hold off;

%%
figure;

plot(splitting_n_mod, means_cauchy);
hold on;

plot(splitting_n_mod, ones(1, n) * a);

xlabel('n');
ylabel('S_n');

legend('Empirical means', 'Shift');

grid on;

hold off;
%%
figure;

plot(sort(cauchy_means), f_means_cauchy, 'r', 'LineWidth', 1.5);
hold on;

plot(sort(cauchy_means), F_cauchy_theor(sort(cauchy_means)), 'b--', 'LineWidth', 1.5);

grid on;

xlabel('x');
ylabel('y');

legend('Empirical distribution function', 'Theoretical distribution function');

axis([-100, 100, 0, 1]);

hold off;