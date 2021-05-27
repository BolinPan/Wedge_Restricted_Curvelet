function [rCurvelets,non_rCurvelets,allCurveletAnglesInfo] = getRestrictedCurvelet(cs_transform, dataType)
% GETRESTRICTEDCURVELET gets restricted Curvelets in all coarse scales
% according to angle theta and specific structure of Curvelet transform
%
%  rCurvelets = getRestrictedCurvelet(theta,cs_transform)
%
%  INPUTS:
%
%   cs_transform - a struct containing parameters:
%         'imageSize'      - image size
%         'nscales'        - number of scales
%         'nangles_coarse' - number of angles in 2nd coarse scale
%         'restrictType'   - select the wedges exactly in the range, i.e.
%                            include those partially in the range
%         'theta' - the maximum angle that wave front impinges on the detector in
%                   radian (0, pi/2)
%
%  OUTPUTS:
%   rCurvelets - restricted Curvelets in 2nd coarse scale
%   allCurveletAngles - a structure containing angle (0, pi/2) orientation 
%                       in image domain 
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke


% assign parameters
nangles_coarse = cs_transform.nangles_coarse;
nscales   = cs_transform.nscales;
imageSize = cs_transform.imageSize;
thetaD = (cs_transform.theta/pi*180); % in degree

% compute the Curvelet mask 
CMask = fdct_wrapping(zeros(imageSize),1,1,nscales,nangles_coarse);

% get parameters of Cp0_Mask
[X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(CMask,imageSize(1),imageSize(2));

% initialization
allCurveletAnglesInfo = [];
rCurvelets = [];

% loop through all angles in each scales (from 2nd scale)
for s = 2:nscales
    for w = 1:length(CMask{s})
        % compute angle (mod pi/2) in frequency domain
        CurveletAngleFrequency = (-1)*atan(F_rows{s}{w}/F_cols{s}{w});

        % compute angle in image domain in degree
        CurveletAngleImage = (pi/2 - atan(tan(CurveletAngleFrequency)/(imageSize(1)/imageSize(2))))/pi*180;
        %CurveletAngleImage = (pi/2 + atan(tan(-CurveletAngleFrequency))/(imageSize(1)/imageSize(2)))/pi*180;

        % put image angles in domain (0,pi/2)
        if CurveletAngleImage < 90
            CurveletAngleImageT = CurveletAngleImage;
        else
            CurveletAngleImageT = 180 - CurveletAngleImage;
        end
        
        % wedges indices
        allCurveletAnglesInfo{s}(w,1) = w;
        
        % get all Curvelet angles in structure
        allCurveletAnglesInfo{s}(w,2) = CurveletAngleImageT;
        
        % indicate restricted wedge
        switch dataType
            case 'data'
                if (CurveletAngleImageT < thetaD)
                    allCurveletAnglesInfo{s}(w,3) = 0; % removed
                else
                    allCurveletAnglesInfo{s}(w,3) = 1; % not removed
                end
            case 'p0'
                if (CurveletAngleImageT > thetaD)
                    allCurveletAnglesInfo{s}(w,3) = 0; % removed
                else
                    allCurveletAnglesInfo{s}(w,3) = 1; % not removed
                end
        end
    end

    % restricted wedges
    rCurvelets_tmp = find(allCurveletAnglesInfo{s}(:,3) == 0);
    rCurvelets_tmpC = find(allCurveletAnglesInfo{s}(:,3) == 1);
    rCurvelets{s} = rCurvelets_tmp.';
    non_rCurvelets{s} = rCurvelets_tmpC.';
end

