
TRY_RING = false;

N = 10;

if TRY_RING % Look at the N-vertex ring graph

    Adj = circshift(eye(N),1);
    [r,c] = find(Adj);
    J = [r,c, -ones(size(r))];
    wG = [r,c];

    Adj = Adj + Adj';

else % Create a 3-regular graph
    [wG, Adj] = randRegGraph(N,3);


    J = [wG, -ones(size(wG,1))];
end


% Find MaxCut to help inform the choice of the Z-type driver

HC = CreateHamC_MaxCut(N, wG);
MaxCut = max(HC);
MaxCut_degen = nnz(HC==MaxCut);

Ixs = find(HC == MaxCut);
Ix = Ixs(randi(length(Ixs)));
aZ = 1+ flip(de2bi(Ix-1,N));

hZ = (-1).^aZ;

figure(1);
plot(graph(Adj));

fprintf('*** N=%d\n', N);
fprintf('    MaxCut = %d, degen = %d\n', MaxCut, MaxCut_degen);
fprintf('    Using hZ with cut value = %d\n', (sum(Adj(:)) - hZ*Adj*hZ')/4);


%% pauli matrices and total angular momentum operators

sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));
Jx = krondist(sx, N)/2;
Jy = krondist(sy, N)/2;
Jz = krondist(sz, N)/2;
Jtot = (Jx^2+Jy^2+Jz^2);

%% Setup QAOA

[QAOAhelperfcn, HamObj, HamC, HamZ, HamB, EvolC, EvolZ, EvolB] = SetupQMCHams(N, J, hZ);


[V, D] = eigs(HamObj);

D = diag(D);

I_GS = D <= D(1)+1e-10;
QMC_degen = nnz(I_GS);

fprintf('Quantum MaxCut value\n')
fprintf('E0 = %0.6f (degen = %d, J^2 = %0.2f)\n', D(1), QMC_degen, V(:,1)'*Jtot*V(:,1));
fprintf('E1 = %0.6f (J^2 = %0.2f)\n',  D(QMC_degen+1), V(:,QMC_degen+1)'*Jtot*V(:,QMC_degen+1))


%% Running and optimizing QAOA

p = 10;

param0 = 2*rand(p,3)-1;

myfun = @(param) QAOAhelperfcn(p, param);

options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
    'TolX',1e-5,'TolFun',1e-5, 'Algorithm', 'quasi-newton','PlotFcn',{@optimplotfval, @optimplotstepsize},...
    'MaxFunEvals', Inf, 'MaxIter', Inf);

[x, fval] = fminunc(myfun, param0, options);

[F1, Fg1, psiQAOA] = myfun(x);

fprintf('QAOA fval = %0.6f @ p=%d\n', fval, p);
fprintf('\t P_GS = %0.4e \n', sum(abs(psiQAOA'*V(:,I_GS)).^2));

%% Setting up AGM ansatz (see arXiv:2003.14394)

zs = {[0;1],[1;0]};

psi0 = 1;
for ind = 1:N
    psi0 = kron(psi0, zs{aZ(ind)});
end


paulis = {sx, sy, sz};

HamYX = sparse(0);

for ind = 1:size(J,1)
    Ia = J(ind,1);
    Ib = J(ind,2);
    sA = paulis{aZ(Ia)};
    sB = paulis{aZ(Ib)};
    HamYX = HamYX + Ham2LTerm(sA, sB, Ia, Ib, N);
end


%% Running and optimizing the AGM ansatz for comparison

yfun = @(t) AGM_ansatz(N, HamObj, HamYX, t, psi0, J, aZ, paulis);

options.PlotFcns = {};
[theta, fval2] = fminunc(yfun, rand*pi, options);

[F2, Fg2, psi_AGM] = AGM_ansatz(N, HamObj, HamYX, theta, psi0, J, aZ, paulis);

fprintf('AGM ansatz fval = %0.6f @ theta = %0.4f\n', fval2, theta);
fprintf('\t P_GS = %0.4e \n', sum(abs(psi_AGM'*V(:,I_GS)).^2));

%%

function [F, Fgrad, psi] = AGM_ansatz(N, HamObj, HamYX, theta, psi0, J, aZ, paulis)
    psi = psi0;
    for ind = 1:size(J,1)
        Ia = J(ind,1);
        Ib = J(ind,2);
        sA = paulis{aZ(Ia)};
        sB = paulis{aZ(Ib)};
        psi = cos(theta) * psi - 1i*sin(theta)*Ham2LTerm(sA, sB, Ia, Ib, N)*psi;
    end
    F = real(psi'*HamObj*psi);
    Fgrad = real(1i*psi'*(HamYX* HamObj -  HamObj*HamYX)*psi);
end
