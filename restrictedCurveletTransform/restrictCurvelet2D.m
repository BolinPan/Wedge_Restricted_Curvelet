function maskC = restrictCurvelet2D(C,rAngle,ID)
%
% RESTRICTCURVELET2D constructs Curvelet mask for 2D data curvelet representation restricting
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
%   maskC  - Curvelet mask for 2D data restrict curvelet representation
%
%
% Copyright (C) 2020 Bolin Pan


% construct all-1 curvelet with structure of C
maskC = constCurvelet(C,true);

% total scales of curvelet
scales = length(maskC); 


% Set maskC to 0 at scales and angles which should be excluded
switch ID
    case {'explicit'}
        %loop through scales
        for s = 2:scales
            if s == 2
                angles = [rAngle, rAngle+(length(maskC{2})/2)];
            elseif s == 3
                rAngleUp = (2*rAngle(1)-1):(2*rAngle(end)+1);
                angles = [rAngleUp, rAngleUp+length(maskC{2})];
            elseif s == 4
                rAngleUp = (2*rAngle(1)-1):(2*rAngle(end)+1);
                angles = [rAngleUp, rAngleUp+length(maskC{2})];
            end
            %loop through angles. Those are infact wedges of tan(angles)
            for w = angles
                maskC{s}{w}(:) = false;  
            end
        end
    case {'implicit'}
        %loop through scales
        for s = 2:scales
            % length of restricted wedges
            L = length(rAngle{s});
            % remove restricted wedges
            for w = 1:L
                maskC{s}{rAngle{s}(w)}(:) = false;  
            end
        end
end
end

