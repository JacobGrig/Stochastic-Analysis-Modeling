% Task 4

%% initials

clc;

x_0 = 1;
gamma = 1;

d = 200;
step = 10000;

n = 100000;

%% calculation

cauchy_vect = cauchy_sensor(x_0, gamma, 1, n);

sorted_cauchy_vect = sort(cauchy_vect);
sorted_cauchy_rest = sorted_cauchy_vect((sorted_cauchy_vect >= -100 + x_0) & (sorted_cauchy_vect <= 100 + x_0));

F_cauchy_emp = (1 : n) ./ n;
F_cauchy_emp_rest = F_cauchy_emp((sorted_cauchy_vect >= -100 + x_0) & (sorted_cauchy_vect <= 100 + x_0));

F_cauchy_theor = @(x) 1 / pi * atan((x - x_0) ./ gamma) + 0.5;

n_splitting = step : step : n;

times_pair = zeros(1, length(n_splitting));
times_model = zeros(1, length(n_splitting));

for i = n_splitting
    if ~(rem(i, 1000))
        disp(i);
    end
    cur_times_pair = zeros(1, d);
    cur_times_model = zeros(1, d);
    for j = 1 : d
        tic;
        [xi_vect, eta_vect] = norm_pair(i);
        cur_times_pair(j) = toc;
        tic;
        norm_vect = norm_model(2 * i);
        cur_times_model(j) = toc;
    end
    times_pair(i / step) = mean(cur_times_pair);
    times_model(i / step) = mean(cur_times_model);
end
%%
norm_vect = norm_model(n);
% norm_vect = zeros(1, n);
% 
% for i = 1 : n
%     norm_vect(i) = norm_model(1);
% end

f_norm_emp = (1 : n) ./ n;

%% visualization

figure;
plot(sorted_cauchy_rest, F_cauchy_emp_rest, 'r');
hold on;

plot(sorted_cauchy_rest, F_cauchy_theor(sorted_cauchy_rest), 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

figure;

plot(sort(norm_vect), f_norm_emp, 'r');
hold on;
splitting_norm = -10 : 0.1 : 10;
plot(splitting_norm, normcdf(splitting_norm), 'b--', 'LineWidth', 1.5);

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function', 'Theoretical distribution function');

grid on;

hold off;

figure;

%norm_super_vect = [norm_vect; 2 * norm_vect; 3 * norm_vect + 1];

normplot(norm_vect);
hold on;

normplot(2 * norm_vect);
normplot(3 * norm_vect + 1);

legend('N(0,1)', 'N(0, sqrt(2))', 'N(1, sqrt(3))');

hold off;


figure;

plot(n_splitting, times_pair, 'b');
hold on;

plot(n_splitting, times_model, 'r');

xlabel('n');
ylabel('t');

legend('Pair method', 'Von Neumann method');

grid on;

hold off;