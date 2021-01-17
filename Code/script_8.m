% Task 7

%% initials

clc;
clear;

h = 5e-2;
n1 = 1e4;
n2 = 1e4;

number_r = 30;
number_phi = 30;

f = @(x,y) x .^ 2 - y .^ 2;

%% calculation

[x_mesh, y_mesh, solution_1] = mk_EMP(h, f, n1);

%%

[r_mesh, phi_mesh, solution_2] = mk_EMP_mod(number_r, number_phi, f, n2);

%% visualization

figure;
surf(x_mesh, y_mesh, solution_1);

xlabel('x');
ylabel('y');
zlabel('z');

sumsqr = x_mesh .^ 2 + y_mesh .^ 2;
bad_indices = sumsqr > 1;

sol_analytic_1 = x_mesh .^ 2 - y_mesh .^ 2;
sol_analytic_1(bad_indices) = NaN;

figure;
surf(x_mesh, y_mesh, sol_analytic_1);

xlabel('x');
ylabel('y');
zlabel('z');

figure;
surf(x_mesh, y_mesh, solution_1 - sol_analytic_1);

xlabel('x');
ylabel('y');
zlabel('z');

figure;
surf(r_mesh .* cos(phi_mesh), r_mesh .* sin(phi_mesh), solution_2);

xlabel('x');
ylabel('y');
zlabel('z');

%hold on;

sol_analytic_2 = (r_mesh .* cos(phi_mesh)) .^ 2 - (r_mesh .* sin(phi_mesh)) .^ 2;

figure;
surf(r_mesh .* cos(phi_mesh), r_mesh .* sin(phi_mesh), sol_analytic_2);

xlabel('x');
ylabel('y');
zlabel('z');

figure;
surf(r_mesh .* cos(phi_mesh), r_mesh .* sin(phi_mesh), solution_2 - sol_analytic_2);

xlabel('x');
ylabel('y');
zlabel('z');