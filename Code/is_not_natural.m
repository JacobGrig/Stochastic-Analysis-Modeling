function res = is_not_natural(number)
res = isnan(number) | imag(number) | (number < 1) | rem(number, 1);
end