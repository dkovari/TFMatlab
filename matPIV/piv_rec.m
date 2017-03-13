function [vx,vy,x,y,vx_nan,vy_nan] = piv_rec(im1,im2,varargin)
% Recursively calculate PIV vectors
% Inputs:
%   im1, im2:   input images (grayscale matricies)
% Parameters:
%   'startW',w  starting width of correlation window (default: 64)
%   'startH',h  starting height of correlation window (default: 64)
%   'endW',w    ending width of correlation window (default: 16)
%   'endH',h    ending height of correlation window (default: 16)
%   'InitScale',s factor by which to rescale initial image during
%                   calculation of image drift
% Output:
%   vx,vy   difference vectors at locations specified by x,y
%   x,y     position of difference vectors (e.g. x=[0.5,16.5,32.5...])
%   vx_nan,vy_nan
%           difference vectors with NaNs at the location where correlation
%           did not return a useable result.  Values of vx,vy at those
%           locations are interpolated from the previous itteration.
%**************************************************************************
% Requirements:
%   piv(...) included with this library
%   imresize() included in the matlab image processing toolbox
%   imxcorr(...) included with this library
%==========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-13: Initial creation (DTK)
%==========================================================================
p=inputParser;
p.CaseSensitive = false;
addParameter(p,'startW',64,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'startH',64,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'endW',16,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'endH',16,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'InitScale',4,@(x) isscalar(x)&&isnumeric(x)&&x>=1);
addParameter(p,'SkipDrift',false,@isscalar);
parse(p,varargin{:});

mindW = p.Results.endW;
mindH = p.Results.endW;

dW=p.Results.startW;
dH = p.Results.startH;

scl = p.Results.InitScale;

fftw('planner','hybrid');

%caclulate whole-size cross-correlation to estimate image drift/average flow
if ~p.Results.SkipDrift
    im1s = imresize(im1,1/scl);
    im2s = imresize(im2,1/scl);
    [dXY,SNR] = imxcorr(im1s,im2s);
    dXY = fix(scl*dXY);
else
    dXY = [0,0];
    SNR = Inf;
end


if ~isempty(dXY)&&SNR>1.5
	x = 0.5:dW/2:size(im1,2)+0.5;
    y = 0.5:dH/2:size(im1,1)+0.5;
    vx = repmat(dXY(1),[numel(y),numel(x)]);
    vy = repmat(dXY(2),[numel(y),numel(x)]);
    
    [vx,vy,x,y,vx_nan,vy_nan] = piv(im1,im2,dW,dH,x,y,vx,vy);
else
    [vx,vy,x,y,vx_nan,vy_nan] = piv(im1,im2,dW,dH);
end

if dW>mindW
    dW = max(mindW,dW/2);
    Wdone = false;
else
    Wdone = true;
end
if dH>mindH
    dH = max(mindH,dH/2);
    Hdone = false;
else
    Hdone = true;
end
while ~Wdone||~Hdone
    %calculate piv with overlap, interpolation
    [vx,vy,x,y,vx_nan,vy_nan] = piv(im1,im2,dW,dH,x,y,vx,vy,'interp');
    %reduce window size
    if dW>mindW
        dW = max(mindW,dW/2);
    else
        Wdone = true;
    end
    if dH>mindH
        dH = max(mindH,dH/2);
    else
        Hdone = true;
    end
end
        