function Q = computeQ(F, A)
% compute modularity
[N, ~] = size(F);
M = sum(sum(A)); % M =  2*m
deg = sum(A, 2);
AD = A - deg * deg' / M;

Ov = sum(F, 2);
for i = 1:N
    F(i, :) = F(i, :) ./ max(Ov(i), 1e-10);
end

tQ = AD .* (F * F');

summ = sum(sum(tQ));

Q = summ / M;

end
