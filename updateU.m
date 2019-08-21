function [Uout] = updateU(A, F, U, Q, IT, alpha)
%%%%%%%%%%%%%%
%% solve ||A - UU^T||^2_F + alpha ||U - FQ||^2_F
%%%%%%%%%%%%%%

i = 1;
fl = 0;
epsilon = 1e-6;

Ucost(1) = f(A, F, U, Q, alpha);

Qpos = (abs(Q) + Q) / 2;
Qneg = (abs(Q) - Q) / 2;
B = 1 / 4;

while (i <= IT) && (fl ~= 1)
    den = (2 * U * (U' * U) + alpha * U + alpha * F * Qneg);
    nom = (2 * A * U + alpha * F * Qpos);
    U = U .* ((nom ./ max(den, 1e-15)).^B);
    
    Ucost(i+1) = f(A, F, U, Q, alpha);
    i = i + 1;
    if (i > 1) && (abs(Ucost(i) - Ucost(i-1)) < epsilon)
        fl = 1;
    end
end
Uout = U;

function c = f(A, F, U, Q, alpha)
tmpA = A - U * U';
tmpF = U - F * Q;
c = sum(diag(tmpA' * tmpA)) + alpha * sum(diag(tmpF' * tmpF));
