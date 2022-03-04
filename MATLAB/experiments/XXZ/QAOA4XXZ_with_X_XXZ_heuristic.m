
N = 12;

maxP = 50;

figind = 11; % index of figure to plot

Delta = -3/4; % gapless = -3/4, gapped = -7/4

p_threshold_for_guess = 6; % needs to be <= 6

nOverlapsToStore = 3; % number of state overlaps to store

ALWAYS_USE_INVERSION_SYM = true;

SAVING_DATA = false;

%% Set up Hamiltonian
sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

HamC = sparse(0);

for ind = 1:N-1
    HamC = HamC + Ham2LTerm(sx, sx, ind, ind+1, N) +...
        Ham2LTerm(sy, sy, ind, ind+1, N) + Delta*Ham2LTerm(sz, sz, ind, ind+1, N);
end

psi0 = ones(2^N,1)/sqrt(2^N);

tic;
if N < 12 && ~ALWAYS_USE_INVERSION_SYM
    [HCsym, Vsym, Dsym] = RestrictToZ2SymSpace(HamC, true);
    fprintf('Finished diagonalizing after %0.2f sec\n', toc);
    
%     Vsub = [speye(2^(N-1)); flipud(speye(2^(N-1)))]/sqrt(2);
%     psi0 = Vsub'*psi0;
    psi0 = ones(2^(N-1),1)/sqrt(2^(N-1));
    
    EvolC = @(psi, gamma) Vsym*(exp(-1i*gamma*Dsym).*(Vsym'*psi));
    EvolB = @(psi, beta) EvolHamB(N, beta, psi, true);
    
    % anonymous function for feeding into optimization (fminunc)
    QAOAhelperfcn = @(p, param) GenQAOAGradSmall(N, p, HCsym, Vsym, Dsym, param, true);
    
    E_GS = Dsym(1); % ground state energy
    E_1E = Dsym(2); % 1st excited state energy
    V_C = Vsym(:, 1:nOverlapsToStore);
else
    % further reduce Hilbert space size using inversion symmetry
    % (but HamB evolution won't work as nicely)
    [HCsym, HBsym, Vsub] = RestrictToZ2SymPlusInversion(N, HamC);
    
    psi0 = Vsub'*psi0;
    
    
    EvolC = @(psi, gamma) expmv(-1i*gamma, HCsym, psi);
    EvolB = @(psi, beta) expmv(-1i*beta, HBsym, psi);
    
    % anonymous function for feeding into optimization (fminunc)
    QAOAhelperfcn = @(p, param) GenQAOAGrad(p, HCsym, HBsym, param, psi0, EvolC, EvolB);
    
    
    % getting lowest eigenstates and eigenvalues
    temp = N^2*speye(length(HCsym)) - HCsym;
    [V_C, Dsym] = eigs(temp, nOverlapsToStore, 'lm');
    
    Dsym = N^2 - diag(Dsym);
    E_GS = Dsym(1); % ground state energy
    E_1E = Dsym(2); % 1st excited state energy
end


%% run QAOA with heuristically guessed initial parameters and educated interpolation strategy

options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
    'TolX',1e-6,'TolFun',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunEvals', Inf, 'MaxIter', Inf);


energyEd = nan(maxP, 1);
paramEd = cell(maxP, 1);
overlapsEd = nan(maxP, nOverlapsToStore);
exitflagEd = nan(maxP, 1);
realTimeEd = nan(maxP, 1);
outputEd = cell(maxP,1);
nFuncEvalEd = nan(maxP, 1);

XXZN8 = load('data/XXZ_n=8_best.mat'); % best known parameters from brute-force search at N=8
% mySymmetry = '';
mySymmetry = 'TR+Z2';

for p = 1:maxP

    if p <= p_threshold_for_guess
        initParam = paramReduce(XXZN8.paramBest{p}, mySymmetry);
    else
        initParam = interpParam(paramReduce(paramEd{p-1}, mySymmetry), p);
    end


    myfun = @(param) QAOAhelperfcn(p, param);

    
    tic;
    [paramEd{p}, Emin, exitflagEd(p), outputEd{p}] = fminunc(myfun, initParam, options);
    realTimeEd(p) = toc;
    nFuncEvalEd(p) = outputEd{p}.funcCount;
    fprintf('p = %d, time = %0.2f sec\n', p, realTimeEd(p));
    
    if p > 1 && Emin > energyEd(p-1)+1e-6 % energy increased relative to previous p
        x = 2;
        while Emin > energyEd(p-1)+1e-6 && p-x >=2
            fprintf('** bad progress, interpolating from p=%d to %d\n', p-x, p);
            tic;
            [param2, Emin2, ef, op] = ...
                fminunc(myfun, interpParam(paramReduce(paramEd{p-x}, mySymmetry), p), options);
            if Emin2 < Emin
                realTimeEd(p) = toc;
                Emin = Emin2;
                paramEd{p} = param2;
                exitflagEd(p) = ef;
                nFuncEvalEd(p) = op.funcCount;
                fprintf('\b + improved, time = %0.2f\n', realTimeEd(p));
            end
            x = x + 1;
        end
    end
    
    
    %%
    energyEd(p) = Emin;
    [~, psiout] = GenQAOA(p, HCsym, paramEd{p}, psi0, EvolC, EvolB);
    
    overlapsEd(p, :) = abs(V_C'*psiout).^2;
    
    figure(figind)
    subplot(2,1,1);
    plot(1:maxP, energyEd,'o-');
    hold on
    plot([1,maxP], E_GS*[1,1], '--r', [1,maxP], E_1E*[1,1], '--k');
    hold off
    xlabel('p');
    ylabel('Energy');
    grid on
    set(gca,'xlim',[1,maxP]);

    subplot(2,1,2);
    plot(1:maxP, overlapsEd, 'o-');
    grid on
    xlabel('p');
    ylabel('eigenstates population');
    legend('Ground','1st Excited', 'location','northwest');
    set(gca,'xlim',[1,maxP]);

    %%
    pause(0.5)
end

%%

figure(figind)
subplot(2,1,1);
plot(1:maxP, energyEd,'o-');
hold on
plot([1,maxP], E_GS*[1,1], '--r', [1,maxP], E_1E*[1,1], '--k');
hold off
xlabel('p');
ylabel('Energy');
grid on
set(gca,'xlim',[1,maxP]);
title(sprintf('N=%d, using interpolation for p>%d',N, p_threshold_for_guess));

subplot(2,1,2);
plot(1:maxP, overlapsEd,'o-');
hold on
% for ind = 1:nOverlapsToStore
%     text(maxP-0.5,mean(overlapsEd(end-1:end, ind)), sprintf('%d',ind-1),'fontsize',16);
% end
hold off, grid on
xlabel('p');
ylabel('eigenstates population');
set(gca,'xlim',[1,maxP]);
set(gca,'ylim',[0,1]);
legend('GS','1E','2E');


%%

if SAVING_DATA
    save(sprintf('data/n=%d_Ed.mat',N), 'maxP','energyEd','overlapsEd','paramEd', 'exitflagEd', 'realTimeEd','nFuncEvalEd', 'E_GS','E_1E');
end
