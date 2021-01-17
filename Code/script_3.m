% Task 3

%% initials

clc;

lambda = 50;
n = 1000;
divider = 0.3;

N = 1000;

alpha = 0.05;

l = 10000;

%% calculation

x_exp_vect = exp_sensor(lambda, n, 1);

sorted_x_exp_vect = sort(x_exp_vect);

F_exp_theory = @(x) 1 - exp(-lambda .* x);

F_exp_emp = (1 : n) / n;

x_exp_vect_rarefied = x_exp_vect(x_exp_vect >= divider);

x_pois1_vect = zeros(1, n);
splitting_n = 1 : n;

for i = splitting_n
    x_pois1_vect(i) = pois_exp_sensor(lambda);
end

pois_theory = @(k) poisspdf(k, lambda);

%%

chi2 = 0;
stud = 0;
fish = 0;

for i = 1 : l
    x_pois2_vect = pois_binom_sensor(lambda, N, n);

    k = max(x_pois2_vect);
    nums = zeros(1, k + 1);
    for j = splitting_n
        cur_elem = x_pois2_vect(j);
        nums(cur_elem + 1) = nums(cur_elem + 1) + 1;
    end

    probs = pois_theory(splitting_n);
    probs_rarefied = pois_theory(0 : k);

    chi_n = n * sum(((nums ./ n - probs_rarefied) .^ 2) ./ probs_rarefied);

    chi_r_theor = chi2inv(1 - alpha, k - 1);

    if (chi_n <= chi_r_theor)
        chi2 = chi2 + 1;
        %disp("Pearson says YES!");
    else
        %disp("Pearson says NO!");
    end
    [xi_vect, eta_vect] = norm_pair(n);
    norm_theor = @(t) exp((-t .^ 2) ./ 2) ./ (sqrt(2 .* pi));

    f_norm_emp = (1 : n) / n;

    xi_mean = mean(xi_vect);
    eta_mean = mean(eta_vect);

    s_1 = sum(((xi_vect - xi_mean).^2) / (n - 1));
    s_2 = sum(((eta_vect - eta_mean).^2) / (n - 1));

    s = (s_1 + s_2) / n;

    t_student = (xi_mean - eta_mean) / sqrt(s);

    t_theory = tinv([alpha/2, 1 - alpha/2], 2 * (n - 2));

    if (t_student <= t_theory(2) && t_student >= t_theory(1))
        stud = stud + 1;
        %disp("Student says YES!");
    else
        %disp("Student says NO!");
    end

    fisher = s_1 / s_2;
    fisher_theor = finv([alpha/2, 1 - alpha/2], n - 1, n - 1);

    if (fisher >= fisher_theor(1) && fisher <= fisher_theor(2))
        fish = fish + 1;
        %disp("Fisher says YES!");
    else
        %disp("Fisher says NO!");
    end
    if (~rem(i,10))
        disp(i);
    end
end

%% visualization

figure;

plot(sorted_x_exp_vect, F_exp_theory(sorted_x_exp_vect), 'r');
hold on;

plot(sorted_x_exp_vect, F_exp_emp, 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Theoretical distribution function', 'Empirical distribution function');

grid on;

hold off;

figure;
hist_1 = histogram(x_exp_vect, 'FaceColor', 'g', 'Normalization', 'probability');
hist_1.BinWidth = 0.005;
hold on;

edges = hist_1.BinEdges;

[x_exp_rarefied, edges_new] = histcounts(x_exp_vect_rarefied - divider, edges, 'Normalization', 'probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;

plot(centers, x_exp_rarefied, 'b', 'LineWidth', 1.5);

xlabel('x');
ylabel('Probability');

legend('Histogram of generated elements', 'Borderline of histogram of rarefied elements, shifted to the left');
grid on;

hold off;

%%

figure;

[x_pois1_hist, edges] = histcounts(x_pois1_vect, 'Normalization','probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;

bar(centers, x_pois1_hist, 2, 'g');
hold on;

bar(splitting_n + 1/4, probs, 0.4, 'b');

%plot(splitting_n, probs, 'b');
%bar(splitting_n, probs, 'b');

xlabel('k');
ylabel('Probability');

legend('Histogram of generated elements', 'Theoretical estimate');
grid on;

axis([20, 80, 0, 0.08]);

hold off;

figure;

[x_pois2_hist, edges] = histcounts(x_pois2_vect, 'Normalization','probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;

bar(centers, x_pois2_hist, 2, 'g');
hold on;

bar(splitting_n + 1/4, probs, 0.4, 'b');

xlabel('k');
ylabel('Probability');

legend('Histogram of generated elements', 'Theoretical estimate');
grid on;

axis([20, 80, 0, 0.08]);

hold off;

figure;

plot(sort(xi_vect), f_norm_emp, 'r');
hold on;
splitting_norm = -10 : 0.1 : 10;
plot(splitting_norm, normcdf(splitting_norm), 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

hold off;

figure;

plot(sort(eta_vect), f_norm_emp, 'r');
hold on;
splitting_norm = -10 : 0.1 : 10;
plot(splitting_norm, normcdf(splitting_norm), 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

hold off;