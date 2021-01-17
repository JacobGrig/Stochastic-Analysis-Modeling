function [r_mesh, phi_mesh, solution] = mk_EMP_mod(number_r, number_phi, f_func, quant)

r_vect = sqrt(linspace(0, 1, number_r));
phi_vect = linspace(0, 2 * pi, number_phi + 1);
phi_vect = phi_vect(1 : (end - 1));

diff_r = [0, diff(r_vect)];
h_phi = phi_vect(2) - phi_vect(1);

[r_mesh, phi_mesh] = meshgrid(r_vect, phi_vect);

solution = nan(number_phi, number_r);

solution(:, end) = f_func(cos(phi_vect), sin(phi_vect));
solution(:, 1) = mean(solution(:, end));

% internal_indices = isnan(solution);
% 
% [init_rows, init_cols] = find(internal_indices);

init_rows = repmat(1 : number_phi, 1, number_r - 2).';
init_cols = repelem(2 : (number_r - 1), number_phi).';

number_of_internals = number_phi * (number_r - 2);

sum_internals = zeros(1, number_of_internals);

for i = 1 : quant
    unhandled_vect = true(1, number_of_internals);
    unhandled = number_of_internals;
    
    rows = init_rows;
    cols = init_cols;
    
    while unhandled
        
        h_i = diff_r(cols);
        h_ip1 = diff_r(cols + 1);
        
        r_i = r_vect(cols);
        
        A = 2 ./ (h_i .* h_ip1) - 1 ./ r_i .* (h_ip1 - h_i) ./ (2 .* h_i .* h_ip1) + 2 ./ (r_i .^ 2) ./ (h_phi .^ 2);
        
        p1 = (2 ./ ((h_i + h_ip1) .* h_ip1) + 1 ./ r_i ./ (2 .* h_ip1)) ./ A;
        p2 = (2 ./ ((h_i + h_ip1) .* h_i) - 1 ./ r_i ./ (2 .* h_i)) ./ A;
        p3 = (1 ./ (r_i .^ 2 * h_phi .^ 2)) ./ A;
        p4 = p3;
        
        steps = count_step([p1; p2; p3; p4]);
%         steps = randi(4, 1, unhandled);
        
        vert_step = zeros(1, unhandled);
        
        vert_step(steps == 3) = 1;
        vert_step(steps == 4) = -1;
        
        horiz_step = zeros(1, unhandled);
        
        horiz_step(steps == 1) = 1;
        horiz_step(steps == 2) = -1;
        
        rows = mod(rows + vert_step.' - 1, number_phi) + 1;
        cols = cols + horiz_step.';
        
        found = (cols == 1) | (cols == number_r);
        
        unhandled = unhandled - sum(found);
        
%         find((cols(found) - 1) * number_phi + rows(found) <= 0)
        
        new_values = solution((cols(found) - 1) * number_phi + rows(found));
        
        inds_to_handled = false(1, number_of_internals);
        inds_to_handled(unhandled_vect) = found;
        
        unhandled_vect(unhandled_vect) = ~found;
        
        sum_internals(inds_to_handled) = sum_internals(inds_to_handled) + new_values.';
        
        rows = rows(~found);
        cols = cols(~found);
    end
    if (~rem(i, 100))
        disp(i);
    end
end

solution((init_cols - 1) * number_phi + init_rows) = sum_internals / quant;
solution = [solution; solution(1, :)];

r_mesh = [r_mesh; r_vect];
phi_mesh = [phi_mesh; 2 * pi * ones(1, number_r)];

end