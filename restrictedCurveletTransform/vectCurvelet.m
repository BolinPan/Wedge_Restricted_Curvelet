function x = vectCurvelet(C, maskC, cSize)
% VECTCURVELET Vectorizes Curvelet coefficients, ordering: scale, angle, (:)
%  x = vectCurvelet(C)
%
% Copyright (C) 2013 Marta M. Betcke

if nargin < 2 || isempty(maskC)
  % If no mask specified, construct all-1 mask
  maskC = constCurvelet(C,1);
end

% acceleration if cSize is known
if exist('cSize','var')
    x = zeros(cSize);
end

i = 0;
for s = 1:length(C) %loop through scales
  for w = 1:length(C{s}) %loop through angles
    lCsw = nnz( maskC{s}{w}(:) );
    x(i + (1:lCsw)) = C{s}{w}( logical( maskC{s}{w}(:) ) );
    i = i + lCsw;
  end
end
  
x = x.'; %return column vector