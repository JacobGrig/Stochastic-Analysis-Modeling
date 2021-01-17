function value = func_f(x)
value = x(1, :) .^ 3 .* sin(1 ./ x(1, :)) + 10 .* x(1, :) .* x(2, :) .^ 4 .* cos(1 ./ x(2, :));
values_temp = x(1, :) .^ 3 .* sin(1 ./ x(1, :));
value(x(2, :) == 0) = values_temp(x(2, :) == 0);
value(x(1, :) == 0) = 0;
end