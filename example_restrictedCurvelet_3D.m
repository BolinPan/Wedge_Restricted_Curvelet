% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
% This script shows the restricted Curvelet transform on 3D ball phantom.
%
%
% Copy right (C) 2021 Bolin Pan

clear all; close all; clc;


%% generate image
Nx = 64;
Ny = 64;
Nz = 64;

% generate 3D ball phantom
radius = Nx/4;
phantom = makeBall(Nx,Ny,Nz,Nx/2,Ny/2,Nz/2,radius);

% display phantom
figure;
imagesc(max(abs(permute(phantom(:,:,:),[2 3 1])),[],3));colorbar; axis off;axis square;
title('$p0$','Interpreter','latex');colorbar;axis image


%% construct wedge restricted Curvelet transform
cs_transform.class         = 'curvelet'; 
cs_transform.nbscales      = 3; % number of scales
cs_transform.nbdstz_coarse = 16; % coarse angles
cs_transform.type          = 'p0R'; % used for initial pressure

% constructs the sparsifying transform Psi - data
cs_transform.imageSize = size(phantom);
disp(['Constructing sparsifying transform: ' cs_transform.class])
% handle to the constructor function of the tranform
cs_transform.constructor = str2func(['@(T) ' cs_transform.class 'CSTransform3D_(T)']);
% construct sparsifying transform
[cs_transform, cs_transform_tag] = cs_transform.constructor(cs_transform);

% construct restricted Curvelet transform
csType = 'explicit';
rCurvelets = [];
rCurvelets{2} = 1:length(cs_transform.maskC{2}); % remove all wedges in 2nd scale (for visualization purpose)
rCurvelets{3} = [];
cs_transform_ = cs_transform;
cs_transform_.maskC = restrictCurvelet3D(cs_transform_.S,rCurvelets,csType); % construct Curvelet mask 4 5 for 32
[cs_transform_, cs_transform_data_tag_] = cs_transform_.constructor(cs_transform_); % Update transform. This needs to be improved to construct in one go. 


%% curvelet transform 
C = cs_transform.unvectPsi(cs_transform.Psi(phantom)); % standard Curvelet
C_ = cs_transform_.unvectPsi(cs_transform_.Psi(phantom)); % restricted Curvelet

% recover ball from curvelet coefficients
p0Recon = cs_transform.iPsi(cs_transform.vectPsi(C));
p0Recon_ = cs_transform_.iPsi(cs_transform_.vectPsi(C_));


%% display and compare
figure;
subplot(2,2,1);imagesc(max(abs(permute(p0Recon(:,:,:),[2 3 1])),[],3));colorbar; axis off;axis square;
title('$\Psi^T \Psi p0$','Interpreter','latex');colorbar;axis image
subplot(2,2,2);imagesc(max(abs(permute(p0Recon_(:,:,:),[2 3 1])),[],3));colorbar; axis off;axis square;
title('$\tilde{\Psi}^\dagger\tilde{\Psi}p0$','Interpreter','latex');colorbar;axis image
subplot(2,2,3);imagesc(max(abs(permute(p0Recon-phantom(:,:,:),[2 3 1])),[],3));colorbar; axis off;axis square;
title('$\Psi^T \Psi p0 - p0$','Interpreter','latex');colorbar;axis image
subplot(2,2,4);imagesc(max(abs(permute(p0Recon-phantom(:,:,:),[2 3 1])),[],3));colorbar; axis off;axis square;
title('$\tilde{\Psi}^\dagger\tilde{\Psi}p0 - p0$','Interpreter','latex');colorbar;axis image




