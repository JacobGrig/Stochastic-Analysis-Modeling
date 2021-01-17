% Task 11

%% initials

clear;
clc;

lambda = 0.15;
T1 = 10000;
N1 = 1000;

T2 = 100;
lambda_0 = 10;

x_m = 1;
k = 2;
W_0 = 300;
c = 0.2;
lambda_ins = 0.1;

N3 = 1000;
T3 = 10000;


%% calculation

[t_QS, requests, Qs, ns_QS] = QS(lambda, N1, T1);

t_QS_mod = QS_mod(lambda_0, T2);
t_QS_mod_mod = QS_mod_mod(lambda_0, T2);

[t_insur, W_insur] = insur(x_m, k, W_0, c, lambda_ins, N3, T3);

%% visualization

figure;

stairs(t_QS, ns_QS);

xlabel('t');
ylabel('n');

grid on;

figure;

stairs(t_QS_mod, 1 : length(t_QS_mod));
hold on;

plot(t_QS_mod, lambda_0 .* (1 + cos(t_QS_mod)));

xlabel('t');
ylabel('y');

legend('W(t)', '\lambda(t)=\lambda_0(1+cos(t))');

grid on;

figure;

stairs(t_QS_mod_mod, 1 : length(t_QS_mod_mod));
hold on;

plot(t_QS_mod_mod, lambda_0 .* (1 + cos(t_QS_mod_mod)));

xlabel('t');
ylabel('y');

legend('W(t)', '\lambda(t)=\lambda_0(1+cos(t))');

grid on;

figure;

plot(t_insur, W_insur);

xlim([min(t_insur), max(t_insur)]);

grid on;

% figure;
% 
% plot(0.1 : 0.1 : 1000, lambda_0 .* (0.1 : 0.1 : 1000 + sin(0.1 : 0.1 : 1000)));