function res = is_not_prob(number)
res = isnan(number) | imag(number) | (number < 0) | (number > 1);
end