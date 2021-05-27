function constC = constCurvelet(C,const)
% CONSTCURVELET Curvelets of size of C with all coefficients equal to const
%
% x = constCurvelet(C)
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke

constC = C;

% Construct a constant coefficients curvelet with structure of C
for s = 1:length(C) % loop through scales
    for w = 1:length(C{s}) % loop through angles. Those are infact wedges of tan(angles)
        constC{s}{w}(:) = const;
    end
end
