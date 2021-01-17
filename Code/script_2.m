% Task 2

%% initials

clc;

n = 1000;
m = 1000;

l = 1000;

epsilon_x = 1e-6;
epsilon_f = 1e-5;

alpha = 0.05;

%% calculation

kolm = 0;
smir1 = 0;
smir2 = 0;

for i = 1 : l

    [x_vect_1, F_vect_1] = cant_sensor(n, epsilon_x);

    [sorted_x_vect_1, sort_ind_1] = sort(x_vect_1);
    sorted_F_vect_1 = F_vect_1(sort_ind_1);

    F_emp_vect_1 = (1 : n) / n;
    F_emp_vect_1_mod = (0 : (n - 1)) / n;

    D_n = max(abs([sorted_F_vect_1 - F_emp_vect_1, sorted_F_vect_1 - F_emp_vect_1_mod]));
    p_value_1 = 1 - kolmcdf(sqrt(n) * D_n, epsilon_f);

    if (p_value_1 > alpha)
        kolm = kolm + 1;
        %disp("Kolmogorov says YES!");
    else
        %disp("Kolmogorov says NO!");
    end

    [x_vect_2, ~] = cant_sensor(m, epsilon_x);
    [sorted_x_vect_2, sort_ind_2] = sort(1 - x_vect_2);

    common_vect = unique(sort([sorted_x_vect_1, sorted_x_vect_2]));
    [sorted_x_vect_1_unique, ind_1] = unique(sorted_x_vect_1, 'last');
    [sorted_x_vect_2_unique, ind_2] = unique(sorted_x_vect_2, 'last');

    F_emp_vect_2 = (1 : m) / m;

    F_emp_new_1 = interp1(sorted_x_vect_1_unique, F_emp_vect_1(ind_1), common_vect, 'previous');
    F_emp_new_2 = interp1(sorted_x_vect_2_unique, F_emp_vect_2(ind_2), common_vect, 'previous');

    % [rows_1, cols_1] = find(sorted_x_vect_1 <= common_vect.');
    % ind_1 = accumarray(rows_1, cols_1, [], @max);
    % 
    % [rows_2, cols_2] = find(sorted_x_vect_2 <= common_vect.');
    % ind_2 = accumarray(rows_2, cols_2, [], @max);
    % 
    % first_ind = max([find(ind_1, 1, 'first'), find(ind_2, 1, 'first')]);
    % 
    % D_nm = max([abs(F_emp_vect_1(ind_1(first_ind : end)) - F_emp_vect_2(ind_2(first_ind : end))), ...
    %     (ind_1(1) == 0) * ((first_ind - 1) / m), (ind_2(1) == 0) * ((first_ind - 1) / n)]);

    D_nm = max(abs(F_emp_new_1 - F_emp_new_2));

    p_value_2 = 1 - kolmcdf(sqrt(n * m / (n + m)) * D_nm, epsilon_f);

    if (p_value_2 > alpha)
        smir1 = smir1 + 1;
        %disp("Smirnov I says YES!");
    else
        %disp("Smirnov I says NO!");
    end

    [x_vect_3, ~] = cant_sensor(m, epsilon_x);
    [sorted_x_vect_3, sort_ind_3] = sort(x_vect_3 / 3);

    sorted_x_vect_1_rarefied = sorted_x_vect_1(sorted_x_vect_1 < (1 / 3));
    n_new = length(sorted_x_vect_1_rarefied);

    common_vect = unique(sort([sorted_x_vect_1_rarefied, sorted_x_vect_3]));
    [sorted_x_vect_1_unique, ind_1] = unique(sorted_x_vect_1_rarefied, 'last');
    [sorted_x_vect_3_unique, ind_3] = unique(sorted_x_vect_3, 'last');

    F_emp_vect_1_rarefied = (1 : n_new) / n_new;
    F_emp_vect_3 = (1 : m) / m;

    F_emp_new_1 = interp1(sorted_x_vect_1_unique, F_emp_vect_1_rarefied(ind_1), common_vect, 'previous');
    F_emp_new_3 = interp1(sorted_x_vect_3_unique, F_emp_vect_3(ind_3), common_vect, 'previous');

    % common_vect = sort([sorted_x_vect_1_rarefied, sorted_x_vect_3]);
    % 
    % F_emp_vect_1_rarefied = (1 : n_new) / n_new;
    % F_emp_vect_3 = (1 : m) / m;
    % 
    % [rows_1, cols_1] = find(sorted_x_vect_1_rarefied <= common_vect.');
    % ind_1 = accumarray(rows_1, cols_1, [], @max);
    % 
    % [rows_3, cols_3] = find(sorted_x_vect_3 <= common_vect.');
    % ind_3 = accumarray(rows_3, cols_3, [], @max);
    % 
    % first_ind = max([find(ind_1, 1, 'first'), find(ind_3, 1, 'first')]);
    % 
    % D_nm = max([abs(F_emp_vect_1_rarefied(ind_1(first_ind : end)) - F_emp_vect_3(ind_3(first_ind : end))), ...
    %     (ind_1(1) == 0) * ((first_ind - 1) / m), (ind_3(1) == 0) * ((first_ind - 1) / n_new)]);

    D_nm = max(abs(F_emp_new_1 - F_emp_new_3));

    p_value_3 = 1 - kolmcdf(sqrt(n_new * m / (n_new + m)) * D_nm, epsilon_f);

    if (p_value_3 > alpha)
        smir2 = smir2 + 1;
        %disp("Smirnov II says YES!");
    else
        %disp("Smirnov II says NO!");
    end
    if (~rem(i, 100))
        disp(i) 
    end
end
    
splitting_n = 1 : n;

% means = zeros(length(splitting_n), 1);
% vars = zeros(length(splitting_n), 1);

[x_vect_cur, ~] = cant_sensor(length(splitting_n), epsilon_x);
means = cumsum(x_vect_cur) ./ (splitting_n);
% vars = (x_vect_cur - means) .^ 2 ./ (1 : length(splitting_n));

vars = cumsum(x_vect_cur .^ 2) ./ splitting_n - means .^ 2;

% for i = 1 : length(splitting_n)
%     %means(i) = mean(x_vect_cur);
%     vars(i) = var(x_vect_cur(1 : i));
% end
%% visualization

figure;

plot(sorted_x_vect_1, sorted_F_vect_1, 'r');
hold on;

plot(sorted_x_vect_1, F_emp_vect_1, 'b--');

xlabel('x');
ylabel('F(x)');

legend('Theoretical distribution function', 'Empirical distribution function');

grid on;

hold off;

figure;

plot(sorted_x_vect_1, F_emp_vect_1, 'r');
hold on;

plot(sorted_x_vect_2, F_emp_vect_2, 'b--');

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function of X', 'Empirical distribution function of (1-X)');

grid on;

hold off;

figure;

plot(sorted_x_vect_1_rarefied, F_emp_vect_1_rarefied, 'r');
hold on;

plot(sorted_x_vect_3, F_emp_vect_3, 'b--');

xlabel('x');
ylabel('F(x)');

legend('Empirical distribution function of X with condition "X<1/3"', 'Empirical distribution function of X/3');

grid on;

hold off;

figure;

semilogx(splitting_n, means, 'r');
hold on;

semilogx(splitting_n, ones(1, length(splitting_n)) / 2, 'b--', 'LineWidth', 2);

xlabel('n');
ylabel('E\xi');

legend('Empirical mean', 'Exact mean');

axis([1, n, -1, 1]);

grid on;

hold off;

figure;

semilogx(splitting_n(2 : end), vars(2 : end), 'r');
hold on;

semilogx(splitting_n(2 : end), ones(1, length(splitting_n) - 1) / 8, 'b--', 'LineWidth', 2);

xlabel('n');
ylabel('D\xi');

legend('Empirical var', 'Exact var');

axis([2, n, -1, 1]);

grid on;

hold off;

