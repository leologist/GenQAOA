N = 10;

maxP = 40;
figind = 21;
p_threshold_for_guess = 9;

SAVING_DATA = false;

%% set up
nOverlapsToStore = 3;

BFseeds = load('data/XXZ_X+aZ+ZZ_n=10_best.mat'); % best known parameters from brute-force search at N=10

[QAOAhelperfcn, HamObj] = SetupXXZHams(N, Inf, true, -3/4);

[V, D] = eig(full(HamObj));
D = diag(D);

V_low = V(:,1:3);
V_low(:, 2:3) = V_low(:, 2:3) *hadamard(2)/sqrt(2);
D_low = D(1:3);

E_GS = D(1);
E_1E = D(2);

myparamReduce = @(param) paramReduceXZZZ(param);

options = optimoptions('fminunc','GradObj','on','Hessian','off','Display','off',...
    'TolX',1e-6,'TolFun',1e-6, 'Algorithm', 'quasi-newton',...
    'MaxFunEvals', Inf, 'MaxIter', Inf);
%%

energyEd = nan(maxP, 1);
paramEd = cell(maxP, 1);
overlapsEd = nan(maxP, nOverlapsToStore);
exitflagEd = nan(maxP, 1);
realTimeEd = nan(maxP, 1);
outputEd = cell(maxP,1);
nFuncEvalEd = nan(maxP, 1);

%%
for p = 1:maxP

    if p <= p_threshold_for_guess
        initParam = BFseeds.paramBest{p};
    else
        initParam = interpParamMultiQAOA(myparamReduce(paramEd{p-1}), p);
    end

    fprintf('p = %d, ', p)
    myfun = @(param) QAOAhelperfcn(p, param);

    %% forward
    
    tic;
    [param, Emin, exitflagEd(p), outputEd{p}] = fminunc(myfun, initParam, options);
    realTimeEd(p) = toc;
    nFuncEvalEd(p) = outputEd{p}.funcCount;
    fprintf('p = %d, time = %0.2f sec\n', p, realTimeEd(p));
    
    if p > 1 && Emin > energyEd(p-1)+1e-6 % energy increased relative to previous p
        x = 2;
        while Emin > energyEd(p-1)+1e-6 && p-x >=2
            fprintf('** bad progress (%0.2e), interp from p=%d to %d\n', energyEd(p-1)-Emin, p-x, p);
            tic;
            [param2, Emin2, ef, op] = ...
                fminunc(myfun, interpParamMultiQAOA(myparamReduce(paramEd{p-x}), p), options);
            if Emin2 < Emin
                realTimeEd(p) = toc;
                fprintf('\b + improved by %0.2e, time = %0.2f\n', Emin-Emin2, realTimeEd(p));
                
                Emin = Emin2;
                param = param2;
                exitflagEd(p) = ef;
                nFuncEvalEd(p) = op.funcCount;
            end
            x = x + 1;
        end
    end
    
    %% backwards
    
%     tic;
%     [param, Emin, exitflagEd(p), outputEd{p}] = fminunc(myfun, initParam, options);
%     realTimeEd(p) = toc;
%     nFuncEvalEd(p) = outputEd{p}.funcCount;
%     fprintf('time = %0.2f sec\n', realTimeEd(p));
%     
%     if p < 40 && Emin > energyEd(p)+1e-6 % energy is larger than known optimum
%         x = 2;
%         while Emin > energyEd(p)+1e-6 && x <= 3
%             fprintf('** bad progress (%0.2e), interp from p=%d to %d\n', energyEd(p)-Emin, p+x, p);
%             tic;
%             [param2, Emin2, ef, op] = ...
%                 fminunc(myfun, interpParamMultiQAOA(myparamReduce(paramEd{p+x}), p), options);
%             if Emin2 < Emin
%                 realTimeEd(p) = toc;
%                 fprintf('\b + improved by %0.2e, time = %0.2f\n', Emin-Emin2, realTimeEd(p));
%                 
%                 Emin = Emin2;
%                 param = param2;
%                 exitflagEd(p) = ef;
%                 nFuncEvalEd(p) = op.funcCount;
%             end
%             x = x + 1;
%         end
%     end
    
    
    %%
    if isnan(energyEd(p)) || Emin < energyEd(p)
        energyEd(p) = Emin;
        paramEd{p} = param;
    end
    [~, ~, psiout] = myfun(paramEd{p});
    
    overlapsEd(p, :) = abs(V_low'*psiout).^2;
    %%
    figure(figind)
    subplot(2,1,1);
    plot(1:maxP, energyEd,'o-');
    hold on
    plot([1,maxP], E_GS*[1,1], '--m', [1,maxP], E_1E*[1,1], '--k');
    hold off
    xlabel('p');
    ylabel('Energy');
    grid on
    legend('\langleH_{XXZ}\rangle','E_0','E_1')
    set(gca,'xlim',[1,maxP]);
    set(gca,'xlim',[1,40]);
    title(sprintf('N=%d, \\Delta=%0.2f',N, Delta))

    subplot(2,1,2);
    plot(1:maxP, overlapsEd, 'o-');
    grid on
    xlabel('p');
    ylabel('eigenstates population');
    legend('Ground','1st Excited', '2nd Excited', 'location','best');
    set(gca,'xlim',[1,maxP]);
    set(gca,'xlim',[1,40]);

    %%
    pause(0.5)
end

if SAVING_DATA
    save('data/XXZ_X+aZ+ZZ_n=10_heu.mat','energyEd','paramEd','overlapsEd','exitflagEd','realTimeEd','outputEd','nFuncEvalEd','maxP','p_threshold_for_guess','N')
end