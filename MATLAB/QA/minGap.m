function [sMin, gap] = minGap(HamC,HamB, varargin)
%MINGAP find the minimum gap in adiabatic algorithm
%initial Hamiltonian HamB = \sum_i X_i
%
% H(s=t/T) = f0(s) * HamB + f1(s) * HamC
% C is the problem Hamiltonian
%
% Goal is to reach highest excited state of C
%
% Required Input:
% HamC: problem Hamiltonian
% HamB: HamB = \sum_i X_i
% 
%   Optional Input:
%       fun0 = f0(s), function of interpolation variable s yielding
%              amplitude of B
%       fun1 = f1(s), function of interpolation variable s yielding
%              amplitude of C
%
% Output:
% gap = minimum gap
% sMin = value of s where minimum gap occurs
%
% Requirement:
% 

%% parsing optional parameter input

ip = inputParser;
addOptional(ip, 'fun0',  @(s) 1-s);
addOptional(ip, 'fun1', @(s) s);
addParameter(ip, 'SvarMax',1);
addParameter(ip, 'SvarMin',0);
addParameter(ip, 'flagRobust',0);

parse(ip, varargin{:});

fun0 = ip.Results.fun0;
fun1 = ip.Results.fun1;
SvarMin = ip.Results.SvarMin;
SvarMax = ip.Results.SvarMax;
flagRobust = ip.Results.flagRobust;


%% Find the minimum gap 

options = optimset('Display','off','TolX',10^-15,'TolFun',10^-15, 'MaxFunEvals', inf, 'MaxIter', inf);

[sMin, gap] = fminbnd(@gapFun,SvarMin,SvarMax,options);

% if robust flag is true, run fminbnd again and fmincon to search for another possibly smaller gap 
if flagRobust 
    [sMin2, gap2] = fminbnd(@gapFun,sMin+0.02,SvarMax,options);    
    options = optimoptions(@fmincon,'GradObj','off','Hessian','off','Display','off','TolX',10^-15,'TolFun',10^-15,'Algorithm','interior-point');        
    [sMin3, gap3] = fmincon(@gapFun,0.8,[],[],[],[],sMin+0.02,1,[],options);
    gap = [gap,gap2,gap3];
    sMin = [sMin,sMin2,sMin3];
    [gap,gapInd] = min(gap);
    sMin = sMin(gapInd);
end

% -------------------------- Gap function -------------------------- %
    function gap = gapFun(s)
        
        Ham = fun0(s)*HamB + fun1(s)*HamC;
        Ham = (Ham+Ham')/2;
        En = eigs(Ham, 2, 'la'); % Find the gap from the highest excited state 
        gap = En(1) - En(2);
    end

end











