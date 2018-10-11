function [SpinObs] = MagObs(psi,N)
% This function computes the single site magnetization Sx, Sy, Sz for the state
% N: number of sites, psi should have length 2^N

SpinObs = zeros(N,3); % N site by sx,sy,sz

for ind=1:N
    SpinObs(ind,:) = [psi'*MultSingleSpin(psi,N,ind,1),...
        real(psi'*MultSingleSpin(psi,N,ind,2)),...
        psi'*MultSingleSpin(psi,N,ind,3)];
end

end