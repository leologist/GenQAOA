graph_type = 1; % see below
experiment_type = 2; % see below

figure_index = 11;


if graph_type == 1 % sublattice of square lattice
    blockadeRadius = sqrt(5);
    xy = rand2Dgrid(5, 0.8)/blockadeRadius;
    NNodes = size(xy, 1);
    [wG, Adj] = UnitDiskGraph(xy);
elseif graph_type == 2 % random unit disk graph
    NNodes = 33;
    [wG, Adj, xy] = UnitDiskGraph(NNodes, 5);
end

% different possible experiment_type
EVOL_MIS_OPT_MIS = 1; % evolve with and optimize MIS Hamiltonian
EVOL_RYD_OPT_MIS = 2; % evolve with Rydberg Hamiltonian, but optimize MIS
EVOL_RYD_OPT_RYD = 3; % evolve with and optimize Rydberg Hamiltonian


%% Setting up Hamiltonians, ground state projectors, etc.

[HamC, HamB] = CreateHams_MIS(NNodes, wG);

[HamC_Ryd, HamB_Ryd, independentSets] = CreateHams_MIS_long_range(xy);

if experiment_type == EVOL_MIS_OPT_MIS || experiment_type == EVOL_RYD_OPT_MIS
    MISsize = max(HamC);
    PrGS = (HamC == MISsize); % projector onto ground states
    Pr1E = (HamC == (MISsize-1)); % prrojector onto first excited states
    E0 = MISsize;
else % experiment_type == EVOL_RYD_OPT_RYD
    E0 = max(HamC_Ryd);
    PrGS = (HamC_Ryd == E0);
    E1 = max(HamC_Ryd(~PrGS));
    Pr1E = HamC_Ryd >= E1 & ~PrGS;
end


fprintf('Hilbert space dimension = %d\n', size(HamC,1));
fprintf('Ground state degen = %d\n', nnz(PrGS)); % typically instances with unique ground state are harder
fprintf('max H_Obj = %f\n', E0);

options = optimoptions('fminunc', 'GradObj', 'on', 'Display', 'off', ...
    'TolFun', 1e-6, 'TolX',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunEvals', Inf, 'MaxIter', Inf);

myPs = 3:20;

%% plot the graph and ground state just for fun
figure(figure_index+2)
gplot(Adj, xy, 'o:')
MIS = independentSets(find(PrGS,1,'first'),:);
hold on
plot(xy(MIS,1), xy(MIS, 2), '.k','markersize',40)
hold off, grid on
title('Graph and excited nodes (black) in ground state')

%% simulate QAOA and optimize

ratioOut = nan(length(myPs), 1);
paramOut = cell(length(myPs),1);
overlapsOut = nan(length(myPs), 2);
exitflag = nan(length(myPs), 1);
realTime = nan(length(myPs), 1);
output = cell(length(myPs),1);
nFuncEval = nan(length(myPs), 1);

for indp = 1:length(myPs)
    p = myPs(indp);
    
    switch experiment_type
        case EVOL_MIS_OPT_MIS
            MISQAOA_fun = @(param) MISQAOAGrad(p, HamC, HamB, param);
        case EVOL_RYD_OPT_MIS
            MISQAOA_fun = @(param) MISQAOAGradExt(p, -HamC, HamC_Ryd, HamB_Ryd, param);
        case EVOL_RYD_OPT_RYD
            MISQAOA_fun = @(param) MISQAOAGradExt(p, -HamC_Ryd, HamC_Ryd, HamB_Ryd, param);
        otherwise
            error('Invalid experiment type %d', experiment_type);
    end
    

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
    ratioOut(indp) = -fmax/E0;
    [~, psiout] = MISQAOA(p, HamC, HamB, paramOut{indp});
    
    overlapsOut(indp, :) = [sum(abs(psiout(PrGS)).^2), sum(abs(psiout(Pr1E)).^2)];
    
    %%
    figure(figure_index)
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

%% plot the optimal QAOA parameters

figure(figure_index+1);
clf
for indp = 1:length(myPs)
    
    p = myPs(indp);
    if isempty(paramOut{indp})
        continue
    end
    subplot(3,6,indp)
    plot(linspace(0,1,p-1), paramOut{indp}(1:p-1), 'o--'); %gamma
    hold on
    plot(linspace(0,1,p), paramOut{indp}(p:end), 'o--'); %beta
    hold off, grid on
    title(sprintf('p=%d', p))
end

