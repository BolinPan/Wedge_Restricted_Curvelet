function [cs_transform, cs_transform_tag] = curveletCSTransform2D_(cs_transform)
% CURVELETCSTRANSFORM Constructs CS transform from Curvelet transform: vectorized image <-> vectorized Curvelet coefficients   
%  [cs_transform, cs_transform_tag] = curveletCSTransform(cs_transform)
%
%  INPUTS:
%   cs_tranform        - structure describing the CS transform with fields
%     class          - transform class: 'curvelet'
%     nscales        - number of scales 
%     imageSize      - an array containing the size of the image on which wavelets operate, [N, M]
%     real           - specifies wether the transform is real (1) of complex (0) 
%     finest         - specifies what is used at finest level: 1 - curvelets, 2 - wavelets, 3 - low pass curvelets
%     transImageSize - size of the image for which curvelets will be constructed; only used with finest=3. 
%                      Restriction: transImageSize <= imageSize
%  OUTPUTS:
%   cs_transform       - adds following fields to CS transform structure
%     S              - structure of the curvelet coefficients 
%     Psi            - handle to forward transform
%     iPsi           - handle to inverse transform
%     imagePsi       - handle to function visualizing Psi coefficients as an image
%   cs_transform_tag - unique tag identifying the transform
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke


if cs_transform.real
  real_tag = 'real';
else
  real_tag = 'compl';
end
if cs_transform.finest == 1
  finest_tag = 'curv';
  % Overwrite transImageSize with imageSize field
  cs_transform.transImageSize = cs_transform.imageSize;
elseif cs_transform.finest == 2
  finest_tag = 'wave';
   % Overwrite transImageSize with imageSize field
  cs_transform.transImageSize = cs_transform.imageSize;
elseif  cs_transform.finest == 3
  % Default, transImageSize value is floor(1/2*imageSize);
  if ~isfield(cs_transform, 'transImageSize')
    cs_transform.transImageSize = floor(1/2*cs_transform.imageSize);
  end
  finest_tag = ['lpcurv' num2str(cs_transform.transImageSize(1)) 'x' num2str(cs_transform.transImageSize(1))];
end

if ~isfield(cs_transform,'nangles_coarse')
  % Default curvelet value of number of coarse scale angles
  cs_transform.nangles_coarse = 16;
end

% Tag describing the transform
cs_transform_tag = [cs_transform.class  '_l' num2str(cs_transform.nscales) '_' real_tag '_f' finest_tag];

% 0-image
sensor0 = zeros(cs_transform.imageSize);

% Construct a dummy curvelet of 0-image
S = fdct_wrapping(sensor0, cs_transform.real, cs_transform.finest, cs_transform.nscales, cs_transform.nangles_coarse);
% Structure of the Curvelet 
cs_transform.S = S;

% Construct all-1 mask if no mask
if ~isfield(cs_transform,'maskC')
  % Default curvelet value of number of coarse scale angles
  cs_transform.maskC = constCurvelet(S, true);
end

% Display of curvelet coefficient (scaled to max 1 at each 'scale' or not)
if ~isfield(cs_transform,'dispScaled')
  % Default curvelet coefficients are displayed scaled to max 1 at each scale
  cs_transform.dispScaled = true;
end

% Construct vectorized coefficient size
cs_transform.vcSize = size(vectCurvelet(S, cs_transform.maskC)');

% Curvelet transform: vectorized image -> vectorized Curvelet
cs_transform.Psi = @(x) vectCurvelet(fdct_wrapping(reshape(x,cs_transform.imageSize),...
cs_transform.real, cs_transform.finest, cs_transform.nscales, cs_transform.nangles_coarse),cs_transform.maskC,cs_transform.vcSize);

% inverse Curvelet
cs_transform.iPsi = @(x_psi) ifdct_wrapping(unvectCurvelet(x_psi, cs_transform.S, cs_transform.maskC),...
    cs_transform.real, cs_transform.transImageSize(1),cs_transform.transImageSize(2));

% Vectorize / unvectorize functions
cs_transform.vectPsi = @(C) vectCurvelet(C, cs_transform.maskC,cs_transform.vcSize);
cs_transform.unvectPsi = @(vC) unvectCurvelet(vC, cs_transform.S, cs_transform.maskC);

% visualization
cs_transform.dispcoef = @(C) curveletDispcoef(C,cs_transform.nscales);

