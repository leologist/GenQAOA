function [F, F_grad, psi_out] = runQAOA4driverMixGrad(p, N, HamObj, EvolHams, params, Evolfuns)
%runQAOA4driverMixGrad runs the 4-driver QAOA with additional variational
%              parameters for the D driver and initial product state, whose
%              derivatives are computed using the finite-different method,
%              and the remaining derivatives are computed analytically
%
%   The 4-drvier QAOA implements the evolution:
%      e^(-1i delta D) e^(-1i gamma C) e^(-1i beta B) e^(-1i alpha A) |m>
%   where
%       A = \sum_ij J_ij Z_i Z_j
%       B = \sum_i X_i
%       C = \sum_i Z_i
%       D = \sum_i n_i . \vec\sigma_i (site-dependent Pauli drivers)
%       |m> = parameterized product state over the qubits
%
%   Usage:
%   [F, F_grad, psi_out] = runQAOA4driverMixGrad(p, N, HamObj,...
%                           EvolHams, params, Evolfuns)
%
%   params = a (4p+4N, 1) vector, where
%           first 2N entries are theta_i, phi_i for n vectors, 
%           next 2N are theta_i, phi_i for m vectorss, 
%           last 4p are alpha1, beta1, gamma1, delta1, alpha2, ... etc.
%


fin_dif = sqrt(eps);
F_grad = zeros(4*N+4*p, 1);

%% construct initial state psi0

ms = reshape(params(2*N+1:4*N), [2, N]).'; % each row is theta_i, phi_i
psi0 = getProductState(ms);

%% construct D driver

n_angles = reshape(params(1:2*N), [2, N]).';
n_thetas = n_angles(:, 1);
n_phis = n_angles(:, 2);
n_vecs = [sin(n_thetas).*cos(n_phis), sin(n_thetas).*sin(n_phis), cos(n_thetas)];

sx = sparse([0,1; 1,0]);
sy = sparse([0,-1i; 1i, 0]);
sz = sparse(diag([1,-1]));
HamD = krondist(sx, N, n_vecs(:,1)) + krondist(sy, N, n_vecs(:,2)) + krondist(sz, N, n_vecs(:,3));

EvolD = @(psi, delta) PauliRotations(N, delta, n_vecs, psi);
EvolHams{end+1} = HamD;
Evolfuns{end+1} = EvolD;


%% get derivatives of the driver parameters

[F, f_grad, psi_out] = MultiQAOAGrad(p, HamObj, EvolHams, reshape(params(4*N+1:end), [4, p]).',...
                                psi0, Evolfuns);

F_grad(4*N+1:end) = reshape(f_grad.', [4*p, 1]);

%% get derivatives of thetas and phis
for i = 1:2*N
    new_pars = params;
    new_pars(i) = new_pars(i) + fin_dif;
    new_angles = reshape(new_pars(1:2*N), [2, N]).';
    new_thetas = new_angles(:, 1);
    new_phis = new_angles(:, 2);
    new_vecs = [sin(new_thetas).*cos(new_phis), sin(new_thetas).*sin(new_phis), cos(new_thetas)];

    EvolD_new = @(psi, delta) PauliRotations(N, delta, new_vecs, psi);
    Evolfuns{end} = EvolD_new;

    f_new = MultiQAOA(p, HamObj, reshape(new_pars(4*N+1:end), [4, p]).', psi0, Evolfuns);
    F_grad(i) = (f_new - F)/fin_dif;
end
Evolfuns{end} = EvolD;
for j = 2*N+1:4*N
    new_pars = params;
    new_pars(j) = new_pars(j) + fin_dif;
    ms = reshape(new_pars(2*N+1:4*N), [2, N]).';
    new_psi0 = getProductState(ms);

    f_new = MultiQAOA(p, HamObj, reshape(new_pars(4*N+1:end), [4, p]).', new_psi0, Evolfuns);
    F_grad(j) = (f_new - F)/fin_dif;
end

end