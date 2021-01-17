function [W_new, R_k] = kalman_filter(W_old, sigma, lambda, sigma_wn, delta_t)

len = length(W_old);
W_new = zeros(1, len + 1);
W_new(1) = 0;

Ck = 1;
Nk = sigma_wn ^ 2;

%%%%%%%
ak = exp(-lambda * delta_t);
Mk = sigma ^ 2 * (1 - exp(-2 * lambda*delta_t));
%%%%%%%

S = sigma ^ 2;

R_k = zeros(1, len + 1);
R_k(1) = S;

for i = 1 : len
    W_kp1_k = ak * W_new(i);
    R_kp1_k = ak * R_k(i) * ak + Mk;
    K = R_kp1_k * Ck / (Ck * R_kp1_k * Ck + Nk) * Ck;
    W_new(i + 1) = W_kp1_k + R_kp1_k * Ck / (Ck * R_kp1_k * Ck + Nk) * (W_old(i) - Ck * W_kp1_k);
    R_k(i + 1) = (1 - K) * R_kp1_k;
end

W_new = W_new(2 : end);
R_k = R_k(2 : end);

end