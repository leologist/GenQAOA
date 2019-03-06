
N = 9;

maxP = 50;

figind = 1; % index of figure to plot

Delta = -3/4; % gapless = -3/4, gapped = -7/4

p_threshold_for_guess = 6; % needs to be <= 6

%% Set up Hamiltonian
sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));

HamC = sparse(0);

for ind = 1:N-1
    HamC = HamC + Ham2LTerm(sx, sx, ind, ind+1, N) +...
        Ham2LTerm(sy, sy, ind, ind+1, N) + Delta*Ham2LTerm(sz, sz, ind, ind+1, N);
end

% work with only first half of the Hilbert space due to X^{\otimes N} symmetry
tic;
[Vsym, Dsym, Hsym] = RestrictToSymSpace(HamC, true);
fprintf('Finished diagonalizing after %0.2f sec\n', toc);


E_GS = Dsym(1); % ground state energy
E_1E = Dsym(2); % 1st excited state energy

options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
    'TolX',1e-6,'TolFun',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunEvals', Inf, 'MaxIter', Inf);

%% educated interpolation strategy

nOverlapsToStore = 7; % number of state overlaps to store

energyEd = nan(maxP, 1);
paramEd = cell(maxP, 1);
overlapsEd = nan(maxP, nOverlapsToStore);
exitflagEd = nan(maxP, 1);
realTimeEd = nan(maxP, 1);
outputEd = cell(maxP,1);
nFuncEvalEd = nan(maxP, 1);


%% run QAOA with heuristically guessed initial parameters

XXZN8 = load('data/XXZ_n=8_best.mat');
% mySymmetry = 'TR+Z2';
mySymmetry = '';

for p = 1:maxP

    if p <= p_threshold_for_guess
        initParam = paramReduce(XXZN8.paramBest{p}, mySymmetry);
    else
        initParam = interpParam(paramReduce(paramEd{p-1}, mySymmetry), p);
    end


    myfun = @(param) GenQAOAGrad(N,p, Hsym, Vsym, Dsym, param, true);

    
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
                fprintf('\b + improved\n');
            end
            x = x + 1;
        end
    end
    
    
    %%
    energyEd(p) = Emin;
    [~, psiout] = GenQAOA(N, p, Hsym, Vsym, Dsym, paramEd{p}, true);
    
    temp = abs(Vsym'*psiout).^2;
    overlapsEd(p, :) = temp(1:nOverlapsToStore);
    
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
    legend('Ground','1st Excited','2nd Excited', 'location','northwest');
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
title(sprintf('N=%d, using interpolation for p>%d',N, p_threshold_for_guess))

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

%%

% save(sprintf('data/n=%d_Ed.mat',N), 'maxP','energyEd','overlapsEd','paramEd', 'exitflagEd', 'realTimeEd','nFuncEvalEd');
