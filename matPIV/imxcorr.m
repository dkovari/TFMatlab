function [dXY,SNR,C]=imxcorr(im1,im2)
% Compute cross-correlation of two images, and return the approximate
% displacement calculated by looking at the shift in peak correlation
% signal.
% Inputs:
%   im1 starting image
%   im2 ending image
% Outputs:
%   dXY, shift in image [x,y] from im1 -> im2
%   SNR     ratio of value of largest peak to second largest peak, or average
%           value of C if no other peaks are detected
%   C   The results of the cross-correlation


% For small windows non-fft based xcorr2 is probably faster.
% O(fft) ~= 12Nlog2(2N)+8N+4
% O(full) ~= 2N^2
% The intersection of the two is N~42.5
if any(size(im1)>42)||any(size(im2)>42)
    C = xcorr2_fft(im2,im1);
else
    C = xcorr2(im2,im1);
end

[dXY,SNR]=imfindpeak(C);

if ~isempty(dXY)
    dXY = dXY-[size(im1,2),size(im1,1)];
end
