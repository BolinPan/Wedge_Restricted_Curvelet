function [u,info] = rL1SALSA(f,para)

% RL1SALSA implements the SALSA algorithm in the reweighted L1 algorithm
% to solve E(x) = 1/2 * || A u - f ||_2^2 + J(u)
% by introducing a split of u as  u - v = 0
% and solving  1/2 * || A u - f ||_2^2 + J(v) such that u - v = 0
% by alternating updates of u, v and w.
% M. V. Afonso, J. M. Bioucas-Dias, and M. A. T. Figueiredo,
% "Fast Image Recovery Using Variable Splitting and Constrainted
% Optimization"
%
%  INPUT:
%   f           - the sensor data f
%   para        - a struct containing all optional parameters:
%                 'overRelaxPara' - over-relaxation parameter (see Section 3.4.3. in
%                  Boyd et al, 2011), default: 1.8, i.e., overrelaxation is in use
%                 'stopCriterion' - choose ... to stop if 'relChangeX' or 'maxIter'
%                 'stopTolerance' - stop tolerance (see above)
%                 'maxIter' - maximum number of iteration after which to stop IN ANY CASE
%                 (even if other convergence criteria are not met yet)
%                 'outputFL' - Logical indicating whether output should be displayed
%                 'tau' - parameter for SALSA
%                 'rho' - parameter for SALSA
%                 'weightMode' - mode for iteratively updating weight
%                 'weightTol'  - tolerance for comparing epsilon
%                 'cs_operator.Pis/iPsi' - forward/backward operators
%
%  OUTPUTS:
%   u           - first version of the primary variable
%   info        - information about the iteration
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke

%%% read out parameters (see above)
stopCriterion     = para.stopCriterion;
stopTolerance     = para.stopTolerance;
maxIter           = para.maxIter;
rho               = para.rho;
tau               = para.tau;
outputFL          = para.outputFL;
overRelaxPara     = para.overRelaxPara;
wSolver           = constructCSWeighting(para);

% assign sparsity operators
PhiPsiT = para.csOperator.Psi;
PsiPhiT  = para.csOperator.iPsi;

% least square solver
hiAATprhoI = @(x,rho) 1/(rho+1)*x;
invLS = @(x,rho) 1/rho*( x - PsiPhiT( hiAATprhoI( PhiPsiT( x ), rho ) ) );

% define operators
E   = @(u) u;
%Etr = @(u) u;
F   = @(v) -v;
b   = zeros(size(PsiPhiT(f)));

%%% initialize inner variables
u0 = b;
u  = u0;
v  = b;
w  = E(u) + F(v) - b;
Fv = F(v);
weight = 1;
eps = 1;

%%% initialize info
iter = 0;
relChangeX = [];
stopFL    = false; 
info = [];

if(outputFL)
    disp('starting reweighted L1 SALSA algorithm iterations.')
end


%%%%%%%%%%%%%%%%%%%%% start the ADMM iteration %%%%%%%%%%%%%%%%%%%%%
while(~stopFL && iter < maxIter)
    
    %%% proceed with the iteration
    iter = iter + 1;
       
    %%% update u
    uOld = u;
    u  = invLS(PsiPhiT(f)+rho*(v - w),rho);
    Eu = E(u);   
    
    %%% over-relaxation 
    auxVar = overRelaxPara * Eu - (1-overRelaxPara) * (Fv - b); % use overRelaxPara = 1
    
    %%% update v via denoising (handling complex number)
    v  = softThresh(auxVar - b + w, tau/rho.*weight);
    
    %%% masking the v 
    Fv = F(v);
   
    %%% update w
    w  = w + (auxVar + Fv - b);
    
    %%% relative change
    changeU   = u - uOld;
    relChangeU = norm(changeU(:))/norm(uOld(:));
    relChangeX(iter,:) = relChangeU;
    
    % update weight
    switch wSolver.mode 
        case 'L0'
            weight = wSolver.weight(u);
        case 'L0r'
            epsOld = eps;
            eps    = wSolver.epsUpdate(u,epsOld);
            weight = wSolver.weight(u,eps);
    end
    
    %%% check stop conditions
    switch stopCriterion
        case 'relChangeX'
            stopValue = relChangeU;
            stopFL = stopValue < stopTolerance;
    end
    
    %%% plotting and output
    switch stopCriterion
        case {'relChangeX','maxIter'}
            outputStr =  ['Iteration ' int2str(iter) '/' int2str(maxIter) ...
                '; relative change : '  num2str(relChangeU,'%.2e')];
    end        
    if(outputFL)
        disp(outputStr)
    end

end


% return some information
info.iter = iter;
info.relChange = relChangeX;

end
