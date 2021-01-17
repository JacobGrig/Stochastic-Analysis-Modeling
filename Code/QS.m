function [t, requests, Qs, ns] = QS(lambda, N, T)

t = 0;

while true
    exps = exprnd(1 / lambda, 1, N);
    new_t = t(end) + cumsum(exps);
    last_ind = find(new_t >= T, 1);
    if ~isempty(last_ind)
        t = [t, new_t(1 : (last_ind - 1))];
        break;
    else
        t = [t, new_t];
    end
end

number_of_requests = length(t);

requests = chi2rnd(10, 1, number_of_requests);

Qs = zeros(1, number_of_requests + 1);
Qs(1) = 0;

splitting = 1 : number_of_requests;

for i = splitting
    if (Qs(i) < t(i))
        Qs(i + 1) = t(i) + requests(i);
    else
        Qs(i + 1) = Qs(i) + requests(i); 
    end
end

Qs = Qs(2 : end);

ns = zeros(1, number_of_requests);

for i = splitting
    ns(i) = sum(Qs(splitting(1 : (i - 1))) > t(i));
end

% ns = sum((splitting.' < splitting) & (Qs(splitting.') > t(splitting)));

end