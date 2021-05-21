function maskC= restrictCurvelet3D(C,rCurvelets,csType)
%
% RESTRICTCURVELET2D constructs Curvelet mask for 3D data curvelet representation restricting
%   curvelet angles to be 0 
%  maskC = restrictCurvelet2D(C,rAngle,ID)
%
%
%  INPUTS:
%   C - a dummy curvelet of 0-image
%   rAngle - restricted angles
%   ID - restrict mode
%
%
%  OUTPUTS:
%   maskC  - Curvelet mask for 3D data restrict curvelet representation
%
%
% Copyright (C) 2020 Bolin Pan


% construct all-1 curvelet with structure of C
maskC = constCurvelet(C,true);

% total scales of curvelet
scales = length(maskC); 


% Set maskC to 0 at scales and angles which should be excluded
switch csType
    case 'explicit'
        for s = 2:scales %loop through scales
            if s == 2
                angles = rCurvelets{2};
            elseif s == 3
                angles = rCurvelets{3};
            elseif s == 4
                angles = rCurvelets{4};
            elseif s == 5
                angles = rCurvelets{5};
            end
            %loop through angles. Those are infact wedges of tan(angles)
            for w = angles
                maskC{s}{w}(:) = false;  
            end
        end
end