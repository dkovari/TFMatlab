function XY = pkfnd_subpixel(im,PKS,method)
% Find subpixel pixel peak
% Input:
%   im grayscale image
%   PKS [X,Y] local maxima locations (e.g. use pkfnd or pkfnd2 to get these)
%   method: 'centroid' use 3-point weighted average (default)
%           'parabolic' estimate peak with parabola
%           'gaussian' use gaussian 3-point fitting
%       Note:   All calculation are currently using 3-point fitting. We might
%               want to change this in the future.
% Output:
%   XY [X,Y] locations of the centroids associated with each peak in PKS
% =========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-12: Initial creation (DTK)
%               All calculations based on 3-point methods presented in:
%                   "Particle Image Velocimetry: A Practical Guide" 2nd Ed
%                   by M. Raffel, et al.
%                   Chapter 5, Table 5.1 pp. 160
%==========================================================================

%% Check method specified
if nargin<3
    method = 'centroid';
end

I = PKS(:,1);
J = PKS(:,2);

X = I;
Y = J;
% remove peaks on the edge
[h,w]=size(im);

ok = find(I~=1&I~=w&J~=1&J~=h);

I = I(ok);
J = J(ok);

%get index from subscript
JI = sub2ind([h,w],J,I);%I+(J-1)*h;%sub2ind([h,w],J,I);

JIn = sub2ind([h,w],J,I-1);%I-1+(J-1)*h;%
JIp = sub2ind([h,w],J,I+1);%I+1+(J-1)*h;%

JnI = sub2ind([h,w],J-1,I);%I+(J-2)*h;%
JpI = sub2ind([h,w],J+1,I);%I+(J)*h;%


switch lower(method)
    case 'centroid'
        X(ok) = ( (I-1).*im(JIn) + I.*im(JI) + (I+1).*im(JIp) )./( im(JIn) + im(JI) + im(JIp) );
        Y(ok) = ( (J-1).*im(JnI) + J.*im(JI) + (J+1).*im(JpI) )./( im(JnI) + im(JI) + im(JpI) );
    case 'parabolic'
        X(ok) = I + ( im(JIn) - im(JIp) )./( 2*im(JIn) - 4*im(JI) +2*im(JIp) );
        Y(ok) = J + ( im(JnI) - im(JpI) )./( 2*im(JnI) - 4*im(JI) +2*im(JpI) );
    case 'gaussian'
        %shift image values so they are always greater than zero
        im([JI;JIn;JIp;JnI;JpI]) = im([JI;JIn;JIp;JnI;JpI]) + 1 - min(im([JI;JIn;JIp;JnI;JpI]));

        %bad = find(im(JIn)==0|im(JI)==0|im(JIp)==0);
        X(ok) = I + ( log(im(JIn)) - log(im(JIp)) )./( 2*log(im(JIn)) - 4*log(im(JI)) +2*log(im(JIp)) );
        %X(ok(bad)) = I(bad) + ( im(JIn(bad)) - im(JIp(ok(bad))) )./( 2*im(JIn(ok(bad))) - 4*im(JI(ok(bad))) +2*im(JIp(ok(bad))) );
        
        %bad = find(im(JnI(ok))==0|im(JI(ok))==0|im(JpI(ok))==0);
        Y(ok) = J + ( log(im(JnI)) - log(im(JpI)) )./( 2*log(im(JnI)) - 4*log(im(JI)) +2*log(im(JpI)) );
        %Y(ok(bad)) = J(ok(bad)) + ( im(JnI(ok(bad))) - im(JpI(ok(bad))) )./( 2*im(JnI(ok(bad))) - 4*im(JI(ok(bad))) +2*im(JpI(ok(bad))) );
    otherwise
        error('Invalid method, must be "centroid" "parabolic" or "gaussian"');
end

XY = [X,Y];