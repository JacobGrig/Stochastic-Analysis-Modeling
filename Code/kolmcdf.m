function f_x = kolmcdf(x, epsilon)
f_x = 1;
k = 1;
mult = -1;
while true
    delta = 2 * mult * exp(-2 * k^2 * x^2);
    if abs(delta) < epsilon
        break;
    end
    f_x = f_x + delta;
    k = k + 1;
    mult = -mult;
end
end