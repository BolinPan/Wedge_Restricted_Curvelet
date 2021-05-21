function img = curveletDispcoef(C, scaling, bkg)
% fdct_wrapping_dispcoef - returns an image containing all the curvelet coefficients
%
% Inputs
%     C         Curvelet coefficients 
%     scale     Scale coefficients at each level to have maximum of 1 (default true)
%     bkg       Value to use for background (default 0.5)
% Outputs
%     img       Image containing all the curvelet coefficients. The coefficents are rescaled so that
%       the largest coefficent in each subband has unit norm.
%
% This function is modified the fdct_wrapping_dispcoef.m in the CurveLab
% toolbox
%
% Copyright (C) 2021 Bolin Pan and Marta M. Betcke  


  if nargin < 3
    bkg = 0.5;
  end
  if nargin < 2
    scaling = true;
  end
  
  [m,n] = size(C{end}{1});
  
  %nbscales = floor(log2(min(m,n)))-3;
  nbscales = length(C);
  if size(C{end},2) == 1
    finest = 2; %Wavelets at finest level
  else
    finest = 1; %Curvelets at finest level
  end
  
  img = C{1}{1};  img = img/max(max(abs(img))); %normalize
  for sc=2:nbscales-(finest-1)
    nd = length(C{sc})/4;
    wcnt = 0;
    
    ONE = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
      ONE = [ONE, fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w})];
    end
    wcnt = wcnt+nd;
    
    TWO = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
      TWO = [TWO; fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w})];
    end
    wcnt = wcnt+nd;
    
    THREE = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
      THREE = [fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w}), THREE];
    end
    wcnt = wcnt+nd;
    
    FOUR = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
      FOUR = [fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w}); FOUR];
    end
    wcnt = wcnt+nd;
    
    [p,q] = size(img);
    [a,b] = size(ONE);
    [g,h] = size(TWO);
    m = 2*a+g;    n = 2*h+b; %size of new image
    if scaling
      scale = max(max( max(max(abs(ONE))),max(max(abs(TWO))) ), max(max(max(abs(THREE))), max(max(abs(FOUR))) )); %scaling factor
      if scale == 0, scale = 1; end
    else
      scale = 1;
    end
    new = bkg * ones(m,n);%background value
    new(a+1:a+g,1:h) = FOUR/scale;
    new(a+g+1:2*a+g,h+1:h+b) = THREE/scale;
    new(a+1:a+g,h+b+1:2*h+b) = TWO/scale;
    new(1:a,h+1:h+b) = ONE/scale;%normalize
    
    dx = floor((g-p)/2);    dy = floor((b-q)/2);
    
    new(a+1+dx:a+p+dx,h+1+dy:h+q+dy) = img;
    
    img = new;
  end

function A = fdct_wrapping_dispcoef_expand(u,v,B)
  A = zeros(u,v);
  [p,q] = size(B);
  A(1:p,1:q) = B;
  
  
  
