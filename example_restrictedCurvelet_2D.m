% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
% This script shows the restricted Curvelet transform on 2D ball phantom.
% 
%
% Copy right (C) 2021 Bolin Pan

clear all; close all; clc;


%% generate image
Nx = 256;
Ny = 256;

% generate 2D ball phantom 
radius = Nx/4;
phantom = makeDisc(Nx,Ny,Nx/2,Ny/2,radius);

% display phantom
figure;
imagesc(phantom);title('$p0$','Interpreter','latex');colorbar;axis image


%% construct wedge restricted Curvelet transform
cs_transform.class          = 'curvelet'; %' mwavelet'; 
cs_transform.nscales        = 3;     % number of scales, higher means more feature windows
cs_transform.nangles_coarse = 16;    % coarse angles, depend on the restricted angle 
cs_transform.real           = 1;     % real curvelets (true/false)
cs_transform.finest         = 1;     % choose 1: curvelets (no doubling),  2: wavelets at finest scale, 3: curvelets (with frequency doubling)
cs_transform.imageSize      = [Nx,Ny];

% handle to the constructor function of the tranform
disp(['Constructing sparsifying transform: ' cs_transform.class])
cs_transform.constructor = str2func(['@(T) ' cs_transform.class 'CSTransform2D_(T)']);
% construct sparsifying transform
[cs_transform, cs_transform_tag] = cs_transform.constructor(cs_transform);

% specify the restricted wedges explicitly in the upper half fourier
% tilling (the angle of wedges ordered from top left closewise)
csType = 'explicit';
upperCurvelets = 1:4; 
cs_transform_ = cs_transform;
% define restircted angles on the frist scale (nagnles_coarse = 128) of upper tilling
cs_transform_.maskC = restrictCurvelet2D(cs_transform_.S,upperCurvelets,csType);
% construct restricted Curvelet transform
[cs_transform_, cs_transform_tag_] = cs_transform_.constructor(cs_transform_);


%% compute curvelet coefficients
C = cs_transform.unvectPsi(cs_transform.Psi(phantom)); % standard Curvelet
C_ = cs_transform_.unvectPsi(cs_transform_.Psi(phantom)); % restricted Curvelet

% display and compare
figure;
subplot(1,2,1);imagesc(abs(cs_transform.dispcoef(C)));title('$|\Psi p0|$','Interpreter','latex')
subplot(1,2,2);imagesc(abs(cs_transform_.dispcoef(C_)));title('$|\tilde{\Psi} p0|$','Interpreter','latex')


%% recover ball phantom from Curvelet domain
p0Recon = cs_transform.iPsi(cs_transform.vectPsi(C));
p0Recon_ = cs_transform_.iPsi(cs_transform_.vectPsi(C));

% display and compare
figure;
subplot(2,2,1);imagesc(p0Recon);axis image;title('$\Psi^T \Psi p0$','Interpreter','latex');colorbar
subplot(2,2,2);imagesc(p0Recon_);axis image;title('$\tilde{\Psi}^\dagger\tilde{\Psi}p0$','Interpreter','latex');colorbar
subplot(2,2,3);imagesc(p0Recon-phantom);axis image;title('$\Psi^T \Psi p0 - p0$','Interpreter','latex');colorbar
subplot(2,2,4);imagesc(p0Recon_-phantom);axis image;title('$\tilde{\Psi}^\dagger\tilde{\Psi}p0 - p0$','Interpreter','latex');colorbar



