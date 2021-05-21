% THIS FILE IS AN EXAMPLE ONLY. DO NOT MODIFY THIS FILE TO RUN EXPERIMENTS!
% MAKE YOUR OWN COPY IN WHICH YOU ADJUST PATHS, OPTIONS ETC.
%
% A script to demonstrate data reconstruction on 2D vessel phantom with
% line sensors geometry in PAT
%
%
% Copy right (C) 2021 Bolin Pan
%
% 

clear all; close all; clc;

% Reset rand, randn, randi to default seed
rng('default')

% define path
path = 'Wedge_Restricted_Curvelet/';

% load setting
load([path,'measurements/sensorMaskUp']) % upscaling sensor mask
load([path,'measurements/subSampleSensorMaskUp']) % upscaling subsampled sensor mask

% set the random generator to the state specified in the noise model
rndSeed = 1;
rng(rndSeed);
sigma = 0.01; % noise level
    
    
%% load phantom
load([path,'measurements/p0'])

% display phantom
figure;imagesc(p0);axis image;title('p0');colorbar


%% run the simulation
source.p0 = p0;

% image size
Nx = size(p0,1);
Ny = size(p0,2);
dx = 1e-4; 
dy = 1e-4; 
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% create the time array
medium.density = 1;
medium.sound_speed = 1500;
kgrid.makeTime(medium.sound_speed);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
dataCast = 'double';
PML_size = [20, 20];
input_args = {'DataCast', dataCast,'PMLSize', PML_size, 'PMLInside', false, 'PlotPML', false, 'Smooth', false, 'PlotSim', false};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add noise
noise = randn(size(sensor_data));
noise = bsxfun(@times, noise, sigma);
fNoisy = sensor_data + noise;


%% get subsampled data
subSampleSensorMaskUp.compressionOperator = [];
subSampleSensorMaskUp.compressionOperator.hA = @(x)full(subSampleSensorMaskUp.A*double(x));
subSampleSensorMaskUp.compressionOperator.hAt = @(x)full(subSampleSensorMaskUp.A'*double(x));
fNoisySS = subSampleSensorMaskUp.compressionOperator.hA(fNoisy);

% display full data and subsampled data
figure;
subplot(3,1,1);imagesc(fNoisy);title('g')
subplot(3,1,2);imagesc(subSampleSensorMaskUp.compressionOperator.hAt(fNoisySS));title('\Phi^T b')
subplot(3,1,3);imagesc(subSampleSensorMaskUp.compressionOperator.hAt(fNoisySS)-fNoisy);title('\Phi^T b - g')


%% compute angle beta based on the bow-tie in Fourier domain
% reorder data before FFT
fNoisyRF = fNoisy.';

% mirror the time domain data about t = 0 to allow the cosine transform to
% be computed using an FFT (p_ty)
fNoisyRF = [flipdim(fNoisyRF, 1); fNoisyRF(2:end, :)];

% extract the size of mirrored input data
[Nt, Ny] = size(fNoisyRF);

% computer kgrid
kgrid = kWaveGrid(Nt, kgrid.dt * medium.sound_speed, Ny, dy);

% compute the angle indicates the inhomogeneous part
beta = atan(kgrid.ky_max/kgrid.kx_max);
disp(['angle beta: ' num2str(beta/pi*180) ' degree'])


%% construct curvelet system for PAT data
cs_transform.class          = 'curvelet'; %' mwavelet'; 
cs_transform.nscales        = 4;     % number of scales, higher means more feature windows
cs_transform.nangles_coarse = 152;   % coarse angles, depend on the restricted angle 
%cs_transform.type           = 'db4'; % type of wavelets 'db1' for Haar, 'db4' etc.
cs_transform.real           = 1;     % real curvelets (true/false)
cs_transform.finest         = 1;     % choose 1: curvelets (no doubling),  2: wavelets at finest scale, 3: curvelets (with frequency doubling)
cs_transform.imageSize      = size(fNoisy);
cs_transform.theta          = beta; % angle for coarse scale visible filter
cs_transform.coarseVisible  = false; % standard Curvelet transform

% obtain restricted Curvelets and angle information 
dataType = 'data';
[rCurvelets, non_rCurvelets, allCurveletAnglesInfo] = getRestrictedCurvelet(cs_transform, dataType);

% handle to the constructor function of the tranform
disp(['Constructing sparsifying transform (data): ' cs_transform.class])
cs_transform.constructor = str2func(['@(T) ' cs_transform.class 'CSTransform2D_(T)']);
% construct sparsifying transform
[cs_transform, cs_transform_tag] = cs_transform.constructor(cs_transform);

% update to restricted Curvelet transform
csType = 'implicit';
cs_transform_ = cs_transform;
% define restircted angles on the frist scale (nagnles_coarse = 128) of upper tilling
cs_transform_.maskC = restrictCurvelet2D(cs_transform_.S,rCurvelets,csType);
% construct restricted Curvelet transform
[cs_transform_, cs_transform_tag_] = cs_transform_.constructor(cs_transform_);

% compute the Curvelet coefficients 
fNoisyC  = cs_transform.Psi(fNoisy); % standard
fNoisyC_ = cs_transform_.Psi(fNoisy); % restricted

% compare the Curvelet coefficients
fNoisyR = cs_transform_.iPsi(fNoisyC_);

% compare the error
figure;
subplot(3,1,1);imagesc(fNoisy);axis image;title('g');colorbar
subplot(3,1,2);imagesc(fNoisyR);axis image;title('g_{r}');colorbar
subplot(3,1,3);imagesc(fNoisyR - fNoisy);axis image;title('g_{r} - g');colorbar


%% full data time-reversal reconstruction
source.p0 = 0;

% obtain sensor mask
sensorMaskAux = false(sensorMaskUp.Nxyz); % upscaled sensor mask reconstruction 
sensorMaskAux(1,:) = sensorMaskUp.sensorMask;
sensor.mask = sensorMaskAux;

% input original full data
sensor.time_reversal_boundary_data = fNoisy;

% run time reversal 
p0TR = kspaceFirstOrder2D(sensorMaskUp.kgrid, medium, source, sensor, input_args{:});

% input compressed full data
sensor.time_reversal_boundary_data = fNoisyR;

% run time reversal 
p0TRR = kspaceFirstOrder2D(sensorMaskUp.kgrid, medium, source, sensor, input_args{:});

% non-negativity processing
p0TRpp = p0TR;
p0TRRpp = p0TRR;
p0TRpp(p0TRpp<0) = 0;
p0TRRpp(p0TRRpp<0) = 0;

% display full data reconstructions
figure;
subplot(3,1,1);imagesc(p0TRpp);axis image;title('p0TRpp');colorbar
subplot(3,1,2);imagesc(p0TRRpp);axis image;title('p0TRpp_{r}');colorbar
subplot(3,1,3);imagesc(p0TRRpp-p0TRpp);axis image;title('p0TRpp_{r} - p0TRpp');colorbar


%% subsampled data time-reversal reconstruction
source.p0 = 0;

% obtain sensor mask
sensorMaskAux = false(subSampleSensorMaskUp.Nxyz);
sensorMaskAux(1,:) = subSampleSensorMaskUp.sensorMask;
sensor.mask = sensorMaskAux;

% input subsampled data
sensor.time_reversal_boundary_data = fNoisySS;

% run time reversal 
p0TRSS = kspaceFirstOrder2D(subSampleSensorMaskUp.kgrid, medium, source, sensor, input_args{:});

% input compressed subsampled data
fNoisySSR = subSampleSensorMaskUp.compressionOperator.hA(fNoisyR);
sensor.time_reversal_boundary_data = fNoisySSR;

% run time reversal 
p0TRRSS = kspaceFirstOrder2D(subSampleSensorMaskUp.kgrid, medium, source, sensor, input_args{:});

% non-negativity processing
p0TRSSpp = p0TRSS;
p0TRRSSpp = p0TRRSS;
p0TRSSpp(p0TRSSpp<0) = 0;
p0TRRSSpp(p0TRRSSpp<0) = 0;

% display subsampled data reconstructions
figure;
subplot(3,1,1);imagesc(p0TRSSpp);axis image;title('p0TRSSpp');colorbar
subplot(3,1,2);imagesc(p0TRRSSpp);axis image;title('p0TRRSSpp');colorbar
subplot(3,1,3);imagesc(p0TRRSSpp-p0TRSSpp);axis image;title('p0TRRSSpp - p0TRSSpp');colorbar


%% solve data reconstruction with R-SALSA
% forward and backward curvelet operators
Psi  = @(x) cs_transform_.Psi(x);
PsiT = @(y) cs_transform_.iPsi(y);
    
% forward operator Phi Psi^-1 
PhiPsiT = @(x) subSampleSensorMaskUp.compressionOperator.hA(PsiT(x));
% adjoint operator (Phi Psi^-1)^T = (Psi^-1)^T Phi^T = Psi Phi^T
PsiPhiT = @(y) Psi(subSampleSensorMaskUp.compressionOperator.hAt(y));
    
% construct inverse method
invParaPP                       = [];
invParaPP.stopCriterion         = 'relChangeX';
invParaPP.stopTolerance         = 5e-4;
invParaPP.maxIter               = 100;
invParaPP.csOperator.Psi        = PhiPsiT;
invParaPP.csOperator.iPsi       = PsiPhiT;
invParaPP.outputFL              = true;
invParaPP.weightMode            = 'L0';
invParaPP.weightTol             = ceil(length(fNoisySS(:))/(5*log(length(fNoisy(:)))));
invParaPP.tau                   = 5e-5;
invParaPP.rho                   = 1;
invParaPP.overRelaxPara         = 1; % no over-relaxation

% reweighted l1 SALSA algorithm
[f_DR,fDRInfo] = rL1SALSA(fNoisySS,invParaPP);

% recover the data from Curvelet coefficients
fNoisyDR = PsiT(f_DR);

% compare reconstructed data with full data
figure;
subplot(1,2,1); imagesc(fNoisyDR); title('g_{DR}')
subplot(1,2,2); imagesc(fNoisyDR-fNoisy);title('g_{DR} - g')

% obtain sensor mask
sensorMaskAux = false(sensorMaskUp.Nxyz); % upscaled sensor mask reconstruction 
sensorMaskAux(1,:) = sensorMaskUp.sensorMask;
sensor.mask = sensorMaskAux;

% input original full data
source.p0 = 0;
sensor.time_reversal_boundary_data = fNoisyDR;

% run time reversal 
p0DR = kspaceFirstOrder2D(sensorMaskUp.kgrid, medium, source, sensor, input_args{:});

% post-processing
p0DRpp = p0DR;
p0DRpp(p0DRpp<0) = 0;

% display reconstruction
figure;imagesc(p0DRpp);axis image;title('p0DR');colorbar


