function [x_vect, F_vect] = cant_sensor(quant, epsilon)

N = ceil(log(1 / epsilon) / log(3));

bern_matr = bern_gen(0.5, N, quant);

degrees = -(1 : N).';

x_vect = sum(2 * bern_matr .* (3 .^ degrees));
F_vect = sum(bern_matr .* (2 .^ degrees));

end