% NNodes = 12;
% wG = randRegGraph(NNodes, 3); % create 3-regular graph

load('demo_QAOAMaxCut_instance.mat','NNodes','wG');
maxP = 30;

%% 
HamC = CreateHamC_MaxCut(NNodes, wG);
HamC = HamC(1:end/2);
CutOpt = max(HamC);
PrGS = HamC >= CutOpt - 1e-6;
Cut1E = max(HamC(~PrGS));
Pr1E = (HamC >= Cut1E - 1e-6) & ~PrGS;

%%
options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
    'TolX',1e-6,'TolFun',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf);

ratioOut = nan(maxP, 1);
paramOut = cell(maxP, 1);
overlapsOut = nan(maxP, 2);
exitflag = nan(maxP, 1);
realTime = nan(maxP, 1);
output = cell(maxP,1);
nFuncEval = nan(maxP, 1);

rng(92702991);

maxP_Trunc = maxP; % maximum number of lower frequency components to keep
numRP = 5;

for p = 1:maxP

    if p == 1
        initParam = [0.79,0.4];
    else
        prevParam = paramOut{p-1};
        prevParam = ConvertMaxCutParam(prevParam, 1);
        prevParam = paramReduceWMC(prevParam);
        prevParam = ConvertMaxCutParam(prevParam, -1);

        if p <= maxP_Trunc
            initParam = [prevParam(1:p-1), 0, prevParam(p:end), 0];
        else
            initParam = [prevParam(1:maxP_Trunc), prevParam((maxP_Trunc+1:2*maxP_Trunc))];
        end
    end

    myfun = @(param) IsingQAOAFourierGrad(NNodes, p, HamC, param, true);
        
    tic;
    [paramOut{p}, fmax, exitflag(p), output{p}] = fminunc(myfun, initParam, options);
    nFuncEval(p) = output{p}.funcCount;
    

    %% try random perturbation (RP)
    
    for trial = 1:numRP
        [param, fmax2, ef, op] = fminunc(myfun, initParam + 0.6*normrnd(0,abs(initParam)), options);
        nFuncEval(p) = nFuncEval(p) + op.funcCount;
        if fmax2 < fmax - 1e-6
            fprintf('*** RP try # %d, Improv = %0.2e, time = %0.2f sec\n', trial, fmax-fmax2, toc);
            fmax = fmax2;
            paramOut{p} = param;
            exitflag(p) = ef;
            output{p} = op;
        end
    end
    
    realTime(p) = toc;
    fprintf('p = %d, time = %0.2f sec\n', p, realTime(p));
    %%
    ratioOut(p) = -fmax/CutOpt;
    
    [~, psiout] = IsingQAOA(NNodes, p, HamC, ConvertMaxCutParam(paramOut{p}, 1, p), true);
    
    overlapsOut(p, :) = [sum(abs(psiout(PrGS)).^2), sum(abs(psiout(Pr1E)).^2)];
    
    figure(1)
    subplot(2,1,1);
    semilogy(1:maxP, 1-ratioOut,'o-');
    xlabel('p');
    ylabel('1-ratio');
    grid on

    subplot(2,1,2);
    plot(1:maxP, overlapsOut, 'o-');
    grid on
    xlabel('p');
    ylabel('overlaps with eigenstates');
    legend('Ground','1st Excited','location','northwest');

    %%
    pause(0.5)
end

%%

load('demo_QAOAMaxCut_instance.mat','wG','resultInterp','resultFourier');

figure(2);
subplot(2,2,[1,3]);
set(0,'DefaultAxesFontSize',12);
graphObj = graph(edge2adj(wG));
plot(graphObj, 'LineWidth',graphObj.Edges.Weight*5);
title('Example Graph');

subplot(2,2,2);
plot(1:resultInterp.maxP, 1-resultInterp.ratio,'.-');
hold on
plot(1:resultFourier.maxP, 1-resultFourier.ratio,'o-');
hold off, grid on
legend('Interp Strategy','Fourier Strategy with 5 RP');
xlabel('p');
ylabel('1 - ratio');
set(gca,'yscale','log');
title('Comparison between two optimization strategies');

subplot(2,2,4);
plot(1:resultInterp.maxP, resultInterp.overlaps(:,1),'.-');
hold on
plot(1:resultFourier.maxP, resultFourier.overlaps(:,1),'o-');
hold off, grid on
legend('Interp Strategy','Fourier Strategy with 5 RP','location','northwest');
xlabel('p');
ylabel('ground state overlap');
