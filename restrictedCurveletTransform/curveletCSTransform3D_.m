function [cs_transform, cs_transform_tag] = curveletCSTransform3D_(cs_transform)
% CURVELETCSTRANSFORM_3D constructs CS transform from complex 3D Curvelet transform: 3D image <-> vectorized Curvelet coefficients   
%  [cs_transform, cs_transform_tag] = curveletCSTransform_3D(cs_transform)
%
%   Construct CS transform from 3D Curvelet transform for Data reconstruction 
%   and initial pressure p0 reconstruction
%
%  INPUTS:
%   cs_transform - a structure describing CS transform with fields
%                'class'          - transform class: 'curvelet'
%                'real'           - specifies whether the transform is real (1) or complex (0)
%                'nscales'        - number of scales
%                'finest'         - specifies what is used at finest level: 1 - curvelets, 2 - wavelets, 3 - low pass curvelets 
%                'nbdstz_coarse'  - define wedges number of frequency tiling in 3D ({8} is equivalent to 32 coarse angles in
%                                   level 2 in 2D frequency tiling}
%                'recon'          - reconstruction: {data_recon, p0_recon}
%                'imageSize'      - an array containing the size of the reshaped 3D image on which wavelets operate, [L,M,N]
%                'transImageSize' - size of the image for which curvelets will be constructed; only used with finest=3. 
%                                       Restriction: transImageSize <= imageSize
%
%  OUTPUTS:
%   cs_transform_tag - unique tage identifying the transfrom
%   cs_transform     - adds following fields to CS transform structure
%                    'S'         - structure of the curvelet coefficients 
%                    'maskC'     - Default curvelet value of number of coarse scale angles
%                                  (construct all-1 mask if no mask)
%                    'Psi'       - handle to forward 3D Curvelet transform
%                    'iPsi'      - handle to backward 3D Curvelet transform
%                    'vectPsi'   - handle to vectorize Curvelet coefficients
%                    'unvectPsi' - handle to unvectorize Curvele coefficients
%
%
% Copyright (C) 2021 Bolin Pan and Marta M. Betcke


% **************************** Initialization *************************** %


% Tag describing the transform
cs_transform_tag = [cs_transform.class  '_l' num2str(cs_transform.nbscales)];

% all Curvelets
allcurvelets = 1;

% 3D 0-image 
o_image = zeros(cs_transform.imageSize);

        
m = size(o_image,1);
n = size(o_image,2);
p = size(o_image,3);

% Construct a dummy curvelet of 0-image
S = fdct3d_forward_mex(m,n,p,...
    cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets, double(o_image));

% Structure of the Curvelet 
cs_transform.S = S;

% Construct all-1 mask if no mask
if ~isfield(cs_transform,'maskC')
  % Default curvelet value of number of coarse scale angles
  cs_transform.maskC = constCurvelet(S, true);
end

% Construct vectorized coefficient size
cs_transform.vcSize = size(vectCurvelet(S, cs_transform.maskC)');

% *********************** 3D Curvelet Transform ************************* %

    
% Backward Curvelet transform: vectorized Curvelet coefficients -> images
switch cs_transform.type
        case 'DR'
            % Forward Curvelet transform: reshaped 3D image -> vectorized Curvelet coefficients
%             cs_transform.Psi = @(x) vectCurvelet(fdct3d_forward_mex(m, n, p,...
%                 cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
%                 double(x)),cs_transform.maskC,cs_transform.vcSize);
%             cs_transform.iPsi = @(x_psi) real(double(fdct3d_inverse_mex(...
%                 m, n, p, cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
%                 unvectCurvelet(x_psi, cs_transform.S, cs_transform.maskC))));
            cs_transform.Psi = @(x) vectCurvelet(fdct3d_forward_mex(m, n, p,...
                cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
                double(reshape(x,cs_transform.imageSize))),cs_transform.maskC,cs_transform.vcSize);
            % Backward Curvelet transform for data reconstruction:
            % vectorized data Curvelet coefficients -> [sensors x timesteps] 2D PAT data 
            % obtained from reshaping 3D PAT data [planar sensors x timesteps]
            cs_transform.iPsi = @(x_psi) real(double(reshape(fdct3d_inverse_mex(...
                m, n, p, cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
                unvectCurvelet(x_psi, cs_transform.S, cs_transform.maskC)), m*n, p)));
        case 'p0R'
            % Forward Curvelet transform: reshaped 3D image -> vectorized Curvelet coefficients
            cs_transform.Psi = @(x) vectCurvelet(fdct3d_forward_mex(m, n, p,...
                cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
                double(reshape(x,cs_transform.imageSize))),cs_transform.maskC,cs_transform.vcSize);

            % Backward Curvelet transform for p0 recosntruction:
            % vectorized p0 Curvelet coeffcients -> initial pressure p0
            cs_transform.iPsi = @(x_psi) real(double(reshape(fdct3d_inverse_mex(...
                m, n, p,...
                cs_transform.nbscales, cs_transform.nbdstz_coarse, allcurvelets,...
                unvectCurvelet(x_psi, cs_transform.S, cs_transform.maskC)),...
                m, n, p)));
end

% Vectorize / unvectorize functions
cs_transform.vectPsi = @(C) vectCurvelet(C,cs_transform.maskC,cs_transform.vcSize);
cs_transform.unvectPsi = @(vC) unvectCurvelet(vC, cs_transform.S, cs_transform.maskC);  
        
end