function psi_out = PauliRotations(N, beta, rot_axes, psi_in)
%PauliRotations efficiently computes exp(-1i*beta* \sum_i sigma_i)*psi_in, 
%                taking advantage of the fact that sigma_i^2 = 1
%
%   Here, sigma_i = (X_i, Y_i, Z_i) .* (n_X, n_Y, n_Z)
%   Since term in the sum commutes and we can break them up, making
%   use of exp(-1i beta S) = cos(beta) * I - i sin(beta) *S whenever S^2=1
%
%   Usage: psi_out = PauliRotations(N, beta, axes, psi_in)
%
%   Arguments:
%       N = system size
%       beta = rotation angle
%       axes = an Nx3 matrix specifying rotation axes (n_X, n_Y, n_Z) for
%              each qubit
%       psi_in = the input wavefunction (2^N dimensional complex vector)
%

for ind = 1:N
    psi_in = cos(beta)*psi_in - 1i*sin(beta)* ...
        (rot_axes(ind, 1) * MultSingleSpin(psi_in,N,ind,1) ...
        + rot_axes(ind, 2) * MultSingleSpin(psi_in,N,ind,2) ...
        + rot_axes(ind, 3) * MultSingleSpin(psi_in,N,ind,3));
end

psi_out = psi_in;

end
