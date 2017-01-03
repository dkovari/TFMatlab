function [qxf,qyf,qxn,qyn,Res] = residfilt(qx,qy,resth,varargin)
% [qxf,qyf,qxn,qyn,Res] = residfilt(qx,qy,...)
% Filter vector field <qx,qy> using Westerweel & Scarano outlier detection
% Inputs:
%   qx: x-components of the vector field
%   qy: y-compnenets of the vector field
% Optional Input:
%   resth:  residual threshold value, above which a vector is deemed an
%           outlier and replaced (default 3.0)
% Parameters:
%   'radius',r  number of neighboring vectors to use in determing if a
%               vector is an outlier (default: 1)
%   'error',er  Approximate estimate of data variability (default: 0.1)
%   'method',str:
%       'median' (default) use median vector value to compute residual and
%                replacement vectors
%       'linear' Fit data in neighborhood to a least-squares approximation
%                of a plane (v=ax+by+c).  Use the expected value to
%                computer residuals
% Outputs:
%   qxf, qyf: The filtered vector field
%   qxn, qyn: The vector field with identified outlier points replaced with
%             NaN
%   Res:    The residual value computed at each location.
%
% The detection algorithm is based on the method presented in:
% Westerweel & Scarano. "Universal outlier detection for PIV data",
% Experiments in Fluids (2005) 39: 1096-110 [DOI:10.1007/s00348-005-0016-6]
%=========================================================================
%% The MIT License
% Copyright (c) 2014 Daniel T. Kovari,
% Georgia Institute of Technology, School of Physics
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software")
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%==========================================================================

%% Change Log =============================================================
%   2014-11-19: File creation (DTK)
%==========================================================================
%% Validate inputs
if any(size(qx)~=size(qy))
    error('qx and qy must be same size');
end

if nargin<3
    resth = 3.0;
end

p=inputParser;
p.CaseSensitive = false;
addOptional(p,'resth',3,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'radius',1,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'error',0.1,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'method','median',@ischar);
parse(p,resth,varargin{:});

resth = p.Results.resth;
rad = p.Results.radius;
er = p.Results.error;
method = p.Results.method;


%calculate the neighborhood rows and cols
[rr,cc] = meshgrid([-rad:-1,1:rad],-rad:rad);
rr=[rr;zeros(1,2*rad)];
cc=[cc;[-rad:-1,1:rad]];
rr=reshape(rr,[],1);
cc=reshape(cc,[],1);

%functions to get neighboring index locations
sz=size(qx);
rc2ind = @(R,C) R+(C-1)*sz(1);
nhoodind = @(r,c) rc2ind(r+rr(rr+r>=1&cc+c>=1&rr+r<=sz(1)&cc+c<=sz(2)),...
                         c+cc(rr+r>=1&cc+c>=1&rr+r<=sz(1)&cc+c<=sz(2)));
hoodexist = @(r,c) rr+r>=1&cc+c>=1&rr+r<=sz(1)&cc+c<=sz(2);

if nargout>4
    Res = zeros(size(qx));
    ResOut=true;
else
    ResOut = false;
end

qxf = qx;
qyf = qy;

if nargout >2
    qxn = qx;
    qxn_out=true;
else
    qxn_out = false;
end
if nargout >3
    qyn = qy;
    qyn_out=true;
else
    qyn_out = false;
end

switch lower(method)
    case 'median'
        for thisR=1:sz(1)
            for thisC=1:sz(2)
                HI = nhoodind(thisR,thisC);    %indecies of neighboring pixels
                if numel(HI)<5 %at a corner and there are too few points for median, use mean
                    Umed = [ nanmean(qx(HI)), nanmean(qy(HI)) ];
                    ResHood = sqrt( (qx(HI)-Umed(1)).^2 + (qy(HI)-Umed(2)).^2);
                    Rs = sqrt( (qx(thisR,thisC)-Umed(1)).^2 + (qy(thisR,thisC)-Umed(2)).^2)/(nanmean(ResHood)+er);
                else
                    Umed = [ nanmedian(qx(HI)), nanmedian(qy(HI)) ];
                    ResHood = sqrt( (qx(HI)-Umed(1)).^2 + (qy(HI)-Umed(2)).^2);
                    Rs = sqrt( (qx(thisR,thisC)-Umed(1)).^2 + (qy(thisR,thisC)-Umed(2)).^2)/(nanmedian(ResHood)+er);
                end
                if ResOut
                    Res(thisR,thisC) = Rs;
                end
                
                if Rs>resth
                    if qxn_out
                        qxn(thisR,thisC) = NaN;
                    end
                    if qyn_out
                        qyn(thisR,thisC) = NaN;
                    end
                    qxf(thisR,thisC) = Umed(1);
                    qyf(thisR,thisC) = Umed(2);
                end
            end
        end
    case 'linear'
        A = [ones(size(rr)),rr,cc];
        bx = zeros(numel(rr),1);
        by = zeros(numel(rr),1);
        for thisR=1:sz(1);
            for thisC=1:sz(2)
                %HI = nhoodind(thisR,thisC);
                HE = hoodexist(thisR,thisC);
                HI = rc2ind(thisR+rr(HE),thisC+cc(HE));
                %get outvals
                bx(HE) = qx(HI);
                bx(~HE) = 0;    %elements outside image go to zero
                by(HE) = qy(HI);
                by(~HE) = 0;    %elements outside image go to zero
                x=lscov(A,[bx,by]);   %compute best planar fit for qx=f(x,y),qy=g(x,y)
                %expected vectors
                qehood = [ones(numel(HI),1),rr(HE),cc(HE)]*x;
                ResHood = sqrt( (qx(HI)-qehood(:,1)).^2 + (qy(HI)-qehood(:,2)).^2);
                thisq = [1,0,0]*x;
                Rs = sqrt( (qx(thisR,thisC)-thisq(1)).^2 + (qy(thisR,thisC)-thisq(2)).^2);
                if numel(HI)<5 %at corner use mean instead of med
                    Rs = Rs/(mean(ResHood)+er);
                else
                    Rs = Rs/(median(ResHood)+er);
                end
                if ResOut
                    Res(thisR,thisC) = Rs;
                end
                if Rs>resth
                    if qxn_out
                        qxn(thisR,thisC) = NaN;
                    end
                    if qyn_out
                        qyn(thisR,thisC) = NaN;
                    end
                    qxf(thisR,thisC) = thisq(1);
                    qyf(thisR,thisC) =thisq(2);
                end
            end
        end
    otherwise
        error('unknown method');
end


