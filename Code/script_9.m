% Task 9

%% initials

clc;
clear;

alpha = 0.4;
n = 16;

lambda = 2;
sigma = 3;

T = 10;

d = 50;

beta = 0.05;

%% calculation

[t_vect_viner, W_vect_viner] = viner_traject(alpha, n);

[t_vect_orn, W_vect_orn] = orn_traject(sigma, lambda, n, T);

%% visualization

figure;

plot(t_vect_viner, W_vect_viner, '.', 'MarkerSize', 1);
hold on;

for i = 1 : d / 2
    [t_vect_viner, W_vect_viner] = viner_traject(alpha, n);
    plot(t_vect_viner, W_vect_viner, '.', 'MarkerSize', 1);
end

plot(t_vect_viner, norminv(1 - beta / 2) .* [-1;1] .* sqrt(t_vect_viner), 'k', 'LineWidth', 1.5);

xlabel('t');
ylabel('W_t');

grid on;

figure;

plot(t_vect_orn, W_vect_orn, '.', 'MarkerSize', 1);
hold on;

for i = 1 : d / 2
    [t_vect_orn, W_vect_orn] = orn_traject(sigma, lambda, n, T);
    plot(t_vect_orn, W_vect_orn, '.', 'MarkerSize', 1);
end

plot([t_vect_orn(1), t_vect_orn(end)], -ones(1, 2) .* norminv(1 - beta / 2) * (sigma ^ 2 / lambda), 'k', 'LineWidth', 1.5);
plot([t_vect_orn(1), t_vect_orn(end)], ones(1, 2) .* norminv(1 - beta / 2) * (sigma ^ 2 / lambda), 'k', 'LineWidth', 1.5);

xlabel('t');
ylabel('W_t');

grid on;