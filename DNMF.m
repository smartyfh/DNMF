function [F] = DNMF(A, U, F, Q, S, para)
i = 1;
alpha = para.alpha;
beta = para.beta;
maxiter = para.maxiter;
IT = para.IT;

while i <= maxiter
    U = updateU(A, F, U, Q, 10*IT, alpha);
    fprintf('finish U in round %d\n', i);
    
    F = updateF(F, U, S, Q, IT, alpha, beta);
    
    Q = updateQ(U, F);
    i = i + 1;
end
