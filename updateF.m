function [Fout] = updateF(F, U, S, Q, IT, alpha, beta)
%%%%%%%%%%%%%%
%% solve alpha ||U - FQ||^2_F + beta tr(F^TSF)
%%%%%%%%%%%%%%

i = 1;
fl = 0;
epsilon = 1e-6;

[N, K] = size(F);
Ik = ones(1, K); 
UQ = U * Q';
SS = beta * S + alpha * eye(N);

Fcost(1) = f(F, U, S, Q, alpha, beta);

while (i <= IT) && (fl ~= 1)            
    for ind = 1:N
        tF = F;
        tF(ind,:) = [];
        s = SS(ind,:);
        s(ind) = [];
        q = UQ(ind,:);
        fi = SS(ind, ind) * Ik + 2 * (s * tF - alpha * q);
        F(ind,:) = h(fi);
    end
    Fcost(i+1) = f(F, U, S, Q, alpha, beta);
    i = i + 1;
    if (i > 1) && (abs(Fcost(i) - Fcost(i-1)) < epsilon) 
        fl = 1;
    end
end
Fout = F;

function c = f(F, U, S, Q, alpha, beta)
tmpU = U - F * Q;
c = alpha * sum(diag(tmpU'*tmpU)) + beta * trace(F' * S * F);

function rf = h(fi)
 K = length(fi);
 rf = zeros(1, K);
 rf(fi<=0) = 1;
 rf(fi == min(fi)) = 1;