% Task 1

%% initials

clc;

p = 0.3;
P = 0.95;

n = 50;
m = 100000;

l = 10000;

N = floor(log(1 - P)/log(1 - p)) * 2;
divider = 5;

%% calculation

bern_matr = bern_gen(p, n, m);
p_cur = sum(sum(bern_matr)) / (n * m);

%%

binom_vect = binom_sensor(p, n, m);
binom_theory = @(k) factorial(n) ./ (factorial(k) .* factorial(n - k)) .* p.^k .* (1 - p).^(n - k);
%binom_theory = pdf('Binomial', ...)
splitting_m = 1 : m;
%geom_vect = zeros(1, m);

% for i = splitting_m
%     geom_vect(i) = geom_gen(p, P);
% end

geom_vect = geom_gen(p, P, m);

geom_theory = @(k) (1 - p) .^ k .* p;

geom_vect_rarefied = geom_vect(geom_vect >= divider);

match_vect = match_pen(m);

%%

matches = zeros(1, l);

for i = 1 : l
    match_cur = match_pen(m);
    matches(i) = match_cur(end);
end

f_norm_emp = (1 : l) ./ l;

%% visualization

figure;

[binom_hist, edges] = histcounts(binom_vect, 'Normalization', 'probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;
bar(centers, binom_hist, 1, 'g');
%histogram(binom_vect, 'FaceColor', 'r', 'Normalization', 'probability');
hold on;

splitting_n = 1 : n;
bar(splitting_n, binom_theory(splitting_n), 0.4, 'b');

xlabel('k');
ylabel('Probability');

legend('Histogram of generated elements', 'Theoretical estimate');
grid on;

hold off;

figure;

[geom_hist, edges] = histcounts(geom_vect, 'Normalization', 'probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;
bar(centers, geom_hist, 1, 'g');

hold on;

splitting_geom = 0 : N;
bar(splitting_geom, geom_theory(splitting_geom), 0.4, 'b');

xlabel('k');
ylabel('Probability');

legend('Histogram of generated elements', 'Theoretical estimate');
grid on;

hold off;

figure;
bar(centers, geom_hist, 1, 'g');
hold on;

[geom_hist_rarefied, edges] = histcounts(geom_vect_rarefied - divider, 'Normalization', 'probability');
centers = (edges(1 : (end - 1)) + edges(2 : end)) ./ 2;

bar(centers, geom_hist_rarefied, 0.4, 'b');

xlabel('k');
ylabel('Probability');

legend('Histogram of generated elements', 'Histogram of rarefied elements, shifted to the left');
grid on;

hold off;

figure;

plot(splitting_m, match_vect);

xlabel('i');
ylabel('Y(i)');

grid on;

%%

figure;

plot(sort(matches), f_norm_emp, 'r', 'LineWidth', 1.5);

hold on;
plot(sort(matches), normcdf(sort(matches)), 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

axis([-3, 4, 0, 1]);

hold off;