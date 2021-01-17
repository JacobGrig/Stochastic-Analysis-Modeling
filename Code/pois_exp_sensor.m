function res_elem = pois_exp_sensor(lambda)
res_elem = 0;
sum_cur = 0;
while true
    sum_cur = sum_cur + exp_sensor(lambda, 1, 1);
    if sum_cur > 1
        break
    end
    res_elem = res_elem + 1;
end
end