function epsilon = getEpsilon(C,epsOld,tol,mode)
% GETEPSILON gets the epsilon value used for IRL1 according to tolerance
% and mode
%
% epsilon = getEpsilon(C,epsOld,tol,mode)
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke


%%% compute the xi0 curvelet
i0 = ceil(tol);
x  = sort(abs(C)./max(abs(C)),1,'descend'); % reorder
xi0 = x(i0);

% check whether element is 0
if xi0 == 0
    switch mode
        case 'L0'
            % nothing to do
        case 'L0r'
            xi0 = epsOld;
    end
end

switch mode
    case 'L0'
        epsilon = max(xi0,1e-4);
    case 'L0r'
        epsilon = min(xi0,epsOld);
end
  
