function wSolver = constructCSWeighting(para)
% CONSTRUCTCSWEIGHTING construct the weighting mode for IRL1
%
% wSolver = constructCSWeighting(para)
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke


% determine whether it is dynamic or static
mode   = para.weightMode;
tol    = para.weightTol;

%%% for fixed epsillon
if isfield(para,'epsilon')
    epsilon = para.epsilon;
    switch mode
        case 'L0'
            weight = @(C) 1./(abs(C)+epsilon);
        case 'L0r'
            weight = @(C) 1./sqrt(abs(C).^2+(epsilon).^2);
    end
else
    switch mode
        case 'L0'
            epsUpdate = @(C) getEpsilon(C,[],tol,mode);
            weight = @(C) 1./(abs(C) + epsUpdate(C));
        case 'L0r'
            epsUpdate = @(C,epsOld) getEpsilon(C,epsOld,tol,mode);
            weight = @(C,epsNew) 1./sqrt(abs(C).^2+(epsNew.^2));
    end
end

wSolver.mode = mode;
wSolver.epsUpdate = epsUpdate;
wSolver.weight = weight;
    

end
