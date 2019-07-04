
N = 8;

Delta = -3/4;

%% Set up Hamiltonian
sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

HamObj = sparse(0);

for ind = 1:N-1
    HamObj = HamObj + Ham2LTerm(sx, sx, ind, ind+1, N) + ...
        Ham2LTerm(sy, sy, ind, ind+1, N) + Delta*Ham2LTerm(sz, sz, ind, ind+1, N);
end


HamC = sparse(0);
for ind = 1:N-1
    HamC = HamC + Ham2LTerm(sz, sz, ind, ind+1, N);
end

HamZ = krondist(sz, N);
HamB = krondist(sx, N);
HamZdiag = full(diag(HamZ));
HamCdiag = full(diag(HamC));

EvolB = @(psi, beta) EvolHamB(N, beta, psi, false);
EvolC = @(psi, gamma) exp(-1i*gamma*HamCdiag).*psi;
EvolZ = @(psi, alpha) exp(-1i*alpha*HamZdiag).*psi;

psi0 = ones(2^N,1)/sqrt(2^N);

QAOAhelperfcn = @(p, param) MultiQAOAGrad(p, HamObj, {HamC, HamZ, HamB}, param, psi0, {EvolC, EvolZ, EvolB});

%% Check gradients


p = 3;
initParam = (2*rand(p, 3)-1)*pi/2;
[F, Fgrad] = QAOAhelperfcn(p, initParam);


nFgrad = nan(size(initParam));


temp = false(size(initParam));
dx = 1e-8;

for ind = 1:numel(nFgrad)
    temp2 = temp; temp2(ind) = 1;
    F2 = QAOAhelperfcn(p, initParam + temp2*dx);
%     F3 = QAOAhelperfcn(p, initParam - temp2*dx);
    nFgrad(ind) = (F2 - F)/dx;
%     nFgrad(ind) = (F2 - F3)/dx/2;
    
end


fprintf('Diff between analytic and numerical gradient:\n\t %0.2e\n', norm(nFgrad- Fgrad))

