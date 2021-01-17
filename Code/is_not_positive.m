function res = is_not_positive(number)
res = isnan(number) | imag(number) | (number <= 0);
end