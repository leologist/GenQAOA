load('demo_QAOA_MIS.mat','NNodes','wG');

[HamC, HamB] = CreateHams_MIS(NNodes, wG);

MISsize = max(HamC);
PrGS = (HamC == MISsize);
Pr1E = (HamC == (MISsize-1));

options = optimoptions('fminunc', 'GradObj', 'on', 'Display', 'off', ...
    'TolFun', 1e-6, 'TolX',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunEvals', Inf, 'MaxIter', Inf);


myPs = 3:20;

%%
ratioOut = nan(length(myPs), 1);
paramOut = cell(length(myPs),1);
overlapsOut = nan(length(myPs), 2);
exitflag = nan(length(myPs), 1);
realTime = nan(length(myPs), 1);
output = cell(length(myPs),1);
nFuncEval = nan(length(myPs), 1);

for indp = 1:length(myPs)
    p = myPs(indp);
    
    MISQAOA_fun = @(param) MISQAOAGrad(p, HamC, HamB, param);

    if p == 3
        gammas = [1.75, -1.75]'; betas = [0.2, 1, 0.4]';
    else
        prev_P = myPs(indp-1);
        prev_gammas = paramOut{indp-1}(1:prev_P-1);
        prev_betas = paramOut{indp-1}(prev_P:2*prev_P-1);

        if prev_P > 1
            gammas = interp1(linspace(0,1,prev_P-1), prev_gammas, linspace(0,1,p-1))';
            betas = interp1(linspace(0,1,prev_P), prev_betas, linspace(0,1,p))';
        else
            gammas = ones(p-1, 1)*prev_gammas;
            betas = ones(p, 1)*prev_betas;
        end
    end
    
    tic;
    [paramOut{indp}, fmax, exitflag(indp), output{indp}] = fminunc(MISQAOA_fun, [gammas;betas], options);
    
    realTime(indp) = toc;
    nFuncEval(indp) = output{indp}.funcCount;
    
    fprintf('Done with p=%d, #iter=%d, #fEval=%d, time = %0.2f sec\n', p, ...
        output{indp}.iterations, output{indp}.funcCount, realTime(indp));
    %%
    ratioOut(indp) = -fmax/MISsize;
    [~, psiout] = MISQAOA(p, HamC, HamB, paramOut{indp});
    
    overlapsOut(indp, :) = [sum(abs(psiout(PrGS)).^2), sum(abs(psiout(Pr1E)).^2)];
    
    %%
    figure(11)
    subplot(2,1,1);
    semilogy(myPs, 1-ratioOut,'o-');
    xlabel('p');
    ylabel('1-ratio');
    grid on

    subplot(2,1,2);
    plot(myPs, overlapsOut, 'o-');
    grid on
    xlabel('p');
    ylabel('overlaps with eigenstates');
    legend('Ground','1st Excited','location','northwest');

    %%
    pause(0.5)
end