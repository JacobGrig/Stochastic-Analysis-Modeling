% Task 10

%% initials

clear;
clc;

n = 9;
sigma = 1;
lambda = 0.1;

T = 1;

sigma_wn = 0.5;
gamma_wn = 0.5;

alpha = 0.05;

%% calculation

[t_vect_orn, W_vect_orn] = orn_traject(sigma, lambda, n, T);
white_noise_norm = sigma_wn * randn(1, 2 ^ n + 1);
white_noise_cauchy = cauchy_sensor(0, gamma_wn, 1, 2 ^ n + 1);

W_vect_orn_noised_norm = W_vect_orn + white_noise_norm;
W_vect_orn_noised_cauchy = W_vect_orn + white_noise_cauchy;

delta_t = t_vect_orn(2) - t_vect_orn(1);

[W_vect_new_norm, R_k_norm] = kalman_filter(W_vect_orn_noised_norm, sigma, lambda, sigma_wn, delta_t);
[W_vect_new_cauchy, R_k_cauchy] = kalman_filter(W_vect_orn_noised_cauchy, sigma, lambda, gamma_wn, delta_t);

bord_norm = W_vect_new_norm +[-1; 1] .* norminv(1 - alpha / 2) * sqrt(R_k_norm);
bord_cauchy = W_vect_new_cauchy +[-1; 1] .* norminv(1 - alpha / 2) * sqrt(R_k_cauchy);
% right_bord_norm = W_vect_new_norm + R_k_norm;

%% visualization

figure;

fill([t_vect_orn, rot90(t_vect_orn, 2)], [bord_norm(1, :), rot90(bord_norm(2, :), 2)], 'g');%, t_vect_orn, right_bord_norm, 'g');
hold on;

plot(t_vect_orn, W_vect_orn, 'LineWidth', 1.5);

% plot(t_vect_orn, W_vect_new_norm);

plot(t_vect_orn, W_vect_orn_noised_norm, '.k', 'MarkerSize', 10);

% plot(t_vect_orn, bord_norm, 'g');%, t_vect_orn, right_bord_norm, 'g');

% errorbar(t_vect_orn, W_vect_new_norm, norminv(1 - alpha / 2) * sqrt(R_k_norm));

legend('possible vector area', 'initial vector', 'noised vector');

grid on;

ylim([-1 + min(bord_norm(1, :)), 1 + max(bord_norm(2, :))]);

figure;

fill([t_vect_orn, rot90(t_vect_orn, 2)], [bord_cauchy(1, :), rot90(bord_cauchy(2, :), 2)], 'g');%, t_vect_orn, right_bord_norm, 'g');

hold on;

plot(t_vect_orn, W_vect_orn, 'LineWidth', 1.5);

% plot(t_vect_orn, W_vect_new_norm);

plot(t_vect_orn, W_vect_orn_noised_cauchy, '.k', 'MarkerSize', 10);

% errorbar(t_vect_orn, W_vect_new_norm, norminv(1 - alpha / 2) * sqrt(R_k_norm));

legend('possible vector area', 'initial vector', 'noised vector');

grid on;

ylim([-7 * gamma_wn + mean(bord_cauchy(1, :)), 7 * gamma_wn + mean(bord_cauchy(2, :))]);