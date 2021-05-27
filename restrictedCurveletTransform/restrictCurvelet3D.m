function maskC= restrictCurvelet3D(C,ID)
%
% RESTRICTCURVELET_3D constructs Curvelet mask for 3D data restrict curvelet representation
%  [maskC, maksVC] = restrictCurvelet_3D(C,ID)
%   
%   Construct Curvelet mask for 3D data curvelet representation restricting
%   curvelet angles not to be 0 (due to causality of wave propagation)
%
%   The list of angles are adapted to 3D Curvelet transform of the Palm
%   data in the paper:
%   Bolin Pan, Simon R. Arridge, Felix Lucka, Ben T. Cox, Nam Huynh, Paul C. Beard, Edward Z. Zhang, Marta M. Betcke 
%   "Photoacoustic Reconstruction Using Sparsity in Curvelet Frame: Image versus Data Domain"
%
%  INPUTS:
%   C - a dummy curvelet of 0-image
%   ID - type of curvelet transform
%
%  OUTPUTS:
%   maskC  - Curvelet mask for 3D data restrict curvelet representation
%
%
% Copyright (C) 2021 Bolin Pan & Marta M. Betcke


% construct all-1 curvelet with structure of C
maskC = constCurvelet(C,true);

% total scales of curvelet
scales = length(maskC); 


% Set maskC to 0 at scales and angles which should be excluded
switch ID
    case 'explicit'
        for s = 2:scales %loop through scales
            if s == 2 
                % 1st face
                angle1 = 673:1092;
                % 2nd face
                angle2Sub = 1781:42:3503;
                angle2 = [angle2Sub,   angle2Sub+1, angle2Sub+2, angle2Sub+3, angle2Sub+4,...
                          angle2Sub+5, angle2Sub+6, angle2Sub+7, angle2Sub+8, angle2Sub+9];
                % 4th face
                angle4 = 5965:6384;
                % 5th face
                angle5Sub = 7073:42:8795;
                angle5 = [angle5Sub,   angle5Sub+1, angle5Sub+2, angle5Sub+3, angle5Sub+4,...
                          angle5Sub+5, angle5Sub+6, angle5Sub+7, angle5Sub+8, angle5Sub+9];
                % all restricted wedges
                angles  =  [angle1, angle2, angle4, angle5];
            elseif (s == 3) || (s == 4)
                % 1st face
                angle1 = 2689:4368;
                % 2nd face
                angle2Sub = 7089:84:14061;
                angle2 = [angle2Sub,    angle2Sub+1,  angle2Sub+2,  angle2Sub+3,  angle2Sub+4,...
                          angle2Sub+5,  angle2Sub+6,  angle2Sub+7,  angle2Sub+8,  angle2Sub+9,...
                          angle2Sub+10, angle2Sub+11, angle2Sub+12, angle2Sub+13, angle2Sub+14,...
                          angle2Sub+15, angle2Sub+16, angle2Sub+17, angle2Sub+18, angle2Sub+19];
                % 4th face
                angle4 = 23857:25536;
                % 5th face
                angle5Sub = 28257:84:35229;
                angle5 = [angle5Sub,    angle5Sub+1,  angle5Sub+2,  angle5Sub+3,  angle5Sub+4,...
                          angle5Sub+5,  angle5Sub+6,  angle5Sub+7,  angle5Sub+8,  angle5Sub+9,...
                          angle5Sub+10, angle5Sub+11, angle5Sub+12, angle5Sub+13, angle5Sub+14,...
                          angle5Sub+15, angle5Sub+16, angle5Sub+17, angle5Sub+18, angle5Sub+19];
                % all restricted wedges
                angles  =  [angle1, angle2, angle4, angle5];
            end
            %loop through angles. Those are infact wedges of tan(angles)
            for w = angles
                maskC{s}{w}(:) = false;  
            end
        end
end