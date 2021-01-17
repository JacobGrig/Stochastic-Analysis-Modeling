function res_vect = geom_gen(p, P, quant)

if is_not_prob(p)
    error('Probability as a parameter of distribution should be a real number between zero and one! Be careful!');
end

if is_not_prob(P)
    error('Probability of "success" should be a real number between zero and one! Be careful!');
end

N = ceil(log(1 - P)/log(1 - p));

res_vect = ones(1, quant) * Inf;
counter = 0;

cur_indices = ones(1, quant);
unhandled = quant;

while unhandled
    bern_vect = bern_gen(p, N, unhandled);
    indices_ones = N - floor(log2(bi2de(rot90(bern_vect, 3))));
    res_vect(logical(cur_indices)) = counter * N + indices_ones - 1;
    cur_indices(~isinf(res_vect)) = 0;
    unhandled = sum(cur_indices);
    %res_vect = counter * N + find(bern_vect, 1, 'first') - 1;
    counter = counter + 1;
end

end