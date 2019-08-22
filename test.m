%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% newF represents the hard community membership matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

% % load adjacent matrix
A = csvread('dolphins-edgesMatrix.csv', 1, 1);

% % initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.5;
beta = 0.1;
gamma = 0.001;
k = 5; %% number of communities to detect
IT = 5; %% number of inner iterations
OT = 50; %% number of outer iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(A, 1);
m = sum(sum(A));
threshold = m / ((n-1) * n);
threshold = sqrt(-log(1-threshold));

I = ones(n, 1);
H = eye(n) - I * I' / n;
tmpA = A' * A;
K = diag(tmpA) * I' + I * (diag(tmpA))' - 2 * tmpA;
K = exp(- 0.5 * K);
Khat = H * K * H;
S = H - (Khat + gamma * eye(n)) \ Khat;

i = 1;
para.alpha = alpha;
para.beta = beta;
para.maxiter = OT;
para.IT = IT;

for iter = 1:1
    U = rand(n, k);
%     F = ones(n, k);
    F = rand(n, k);
    Q = rand(k, k);
    Q = ProjTF(Q);
    [newF] = DNMF(A, U, F, Q, S, para);
    modularity = computeQ(newF, A)
end
