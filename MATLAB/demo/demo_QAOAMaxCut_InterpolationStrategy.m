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

figure
plot(graph(edge2adj(wG)));

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

for p = 1:maxP

    if p == 1
        initParam = [0.8,0.35];
    elseif p == 2
        initParam = [paramOut{1}(1)*[1,1], paramOut{1}(2)*[1,1]];
    else
        xp = linspace(0,1,p-1);
        xp1 = linspace(0,1,p);
        initParam = [interp1(xp, paramOut{p-1}(1:p-1), xp1, 'linear'), ...
                     interp1(xp, paramOut{p-1}(p:end), xp1, 'linear')];
    end


    myfun = @(param) IsingQAOAGrad(NNodes, p, HamC, param, true);

    tic;
    [paramOut{p}, fmax, exitflag(p), output{p}] = fminunc(myfun, initParam, options);
    realTime(p) = toc;
    nFuncEval(p) = output{p}.funcCount;
    
    fprintf('p = %d, time = %0.2f sec\n', p, realTime(p));
    %%
    ratioOut(p) = -fmax/CutOpt;
    [~, psiout] = IsingQAOA(NNodes, p, HamC, paramOut{p}, true);
    
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
