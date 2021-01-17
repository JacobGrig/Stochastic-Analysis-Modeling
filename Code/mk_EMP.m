function [x_mesh, y_mesh, solution] = mk_EMP(step, f_func, quant)
x_vect = -1 : step : 1;
y_vect = -1 : step : 1;

len_x = length(x_vect);
len_y = length(y_vect);

[x_mesh, y_mesh] = meshgrid(x_vect, y_vect);

solution = inf(len_x, len_y);

sumsqr = x_mesh .^ 2 + y_mesh .^ 2;

external_indices = sumsqr > 1;

solution(external_indices) = NaN;

border_indices = false(len_x, len_y);

up_indices = sumsqr(1, :) <= 1;
down_indices = sumsqr(end, :) <= 1;

left_indices = sumsqr(:, 1) <= 1;
right_indices = sumsqr(:, end) <= 1;

rest_border = (sumsqr(2 : (end - 1), 2 : (end - 1)) <= 1) & ...
             ((sumsqr(1 : (end - 2), 2 : (end - 1)) > 1) | ...
              (sumsqr(3 : (end - 0), 2 : (end - 1)) > 1) | ...
              (sumsqr(2 : (end - 1), 1 : (end - 2)) > 1) | ...
              (sumsqr(2 : (end - 1), 3 : (end - 0)) > 1));

border_indices(1, up_indices) = true;
border_indices(end, down_indices) = true;
border_indices(left_indices, 1) = true;
border_indices(right_indices, end) = true;
          
border_indices(2 : (end - 1), 2 : (end - 1)) = rest_border;

solution(border_indices) = f_func(x_mesh(border_indices), y_mesh(border_indices));

internal_indices = (~external_indices) & (~border_indices);

[init_rows, init_cols] = find(internal_indices);

number_of_internals = length(init_rows);

sum_internals = zeros(1, number_of_internals);

for i = 1 : quant
    unhandled_vect = true(1, number_of_internals);
    unhandled = number_of_internals;
    
    rows = init_rows;
    cols = init_cols;
    
    while unhandled
        steps = randi(4, 1, unhandled);
        
        vert_step = zeros(1, unhandled);
        
        vert_step(steps == 2) = 1;
        vert_step(steps == 4) = -1;
        
        horiz_step = zeros(1, unhandled);
        
        horiz_step(steps == 1) = -1;
        horiz_step(steps == 3) = 1;
        
        rows = rows + vert_step.';
        cols = cols + horiz_step.';
        
        found = border_indices((cols - 1) * len_x + rows);
        
        unhandled = unhandled - sum(found);
        
        new_values = solution((cols - 1) * len_x + rows);
        
        inds_to_handled = false(1, number_of_internals);
        inds_to_handled(unhandled_vect) = found;
        
        unhandled_vect(unhandled_vect) = ~found;
        
        sum_internals(inds_to_handled) = sum_internals(inds_to_handled) + new_values(found).';
        
        rows = rows(~found);
        cols = cols(~found);
    end
    if (~rem(i, 100))
        disp(i);
    end
end

solution((init_cols - 1) * len_x + init_rows) = sum_internals / quant;

end