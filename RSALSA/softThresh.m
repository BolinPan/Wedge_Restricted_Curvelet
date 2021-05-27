function y = softThresh(x,T)
% SOFTTHRESH Soft thresholding
%
%  softThresh(x,T)
%
%
% Copy right (C) 2021 Bolin Pan & Marta M. Betcke

size_x = size(x);

y = max(abs(x(:)) - T, 0);
y = y./(y+T) .* x(:);

y = reshape(y, size_x);
