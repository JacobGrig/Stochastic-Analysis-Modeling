function res = is_not_real(number)
res = isnan(number) | imag(number);
end