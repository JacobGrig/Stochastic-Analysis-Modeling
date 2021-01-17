function res_vect = norm_model(quant)

if is_not_natural(quant)
    error('An amount of elements should be a natural number! Be careful!');
end

%k = sqrt(2 * pi / exp(1));
gamma = 1;
x_0 = 0;

res_vect = nan(1, quant);

unhandled = quant;

%indices = ones(1, quant);

% p2_p1_func = @(x) sqrt(pi ./ 2) ./ gamma .* (gamma .^ 2 + (x - x_0) .^ 2) .* exp((-x .^ 2) ./ 2);

p2_kp1_func = @(x) exp((1 - x .^ 2) ./ 2) .* (x .^ 2 + 1) ./ 2;

% p1_func = @(x) 1 ./ pi .* gamma ./ ((x - x_0) .^ 2 + gamma .^ 2);
% p2_func = @(x) 1 ./ sqrt(2 * pi) * exp(-(x .^ 2) / 2);

while unhandled
    cauchy_vect = cauchy_sensor(x_0, gamma, 1, unhandled);
    % probs_mod = p2_func(cauchy_vect) ./ (k .* p1_func(cauchy_vect));
    probs = p2_kp1_func(cauchy_vect);
    bern_vect = logical(rand(1, unhandled) < probs);
    
    number_ones = sum(bern_vect);
    res_vect((quant - unhandled + 1) : (quant - unhandled + number_ones)) = cauchy_vect(bern_vect);
    unhandled = unhandled - number_ones;
%     cur_indices = logical(indices);
%     cur_indices(logical(indices)) = bern_vect;
%     
%     res_vect(cur_indices) = cauchy_vect(bern_vect);
%     
%     indices(cur_indices) = 0;
%     unhandled = sum(indices);
end

end