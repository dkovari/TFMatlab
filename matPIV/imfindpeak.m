function [XY,SNR] = imfindpeak(C)
% Finds the highest peak in an image and returns the subpixel location of
% that peak using a 3-point gaussian fitting.
% Input:
%   C input image
% Output:
%   XY = [X,Y] the sub-pixel location of the highest peak
%   SNR     ratio of value of largest peak to second largest peak, or average
%           value of C if no other peaks are detected
%**************************************************************************
% Requirements:
%   pkfnd(...) by Dufresnes et al. (needs DTK modification).
%   pkfnd_subpixel(...) by D. Kovari
%==========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-13: Initial creation (DTK)
%==========================================================================

[pks,ind] = pkfnd(C);

if isempty(ind)
    XY = [];
    SNR = 0;
    return;
end

[mxpkv,mxpk] = max(C(ind));
XY = pks(mxpk,1:2);

ind(mxpk) = [];
mxpkv2 = max(C(ind));
if isempty(mxpkv2)
    mxpkv2 = nanmean(C(:));
end

%calc SNR (need to look this up)
SNR = mxpkv/mxpkv2;
XY = pkfnd_subpixel(C,XY,'gaussian');