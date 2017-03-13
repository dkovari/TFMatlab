function res = bpass(image_array,lnoise,lobject)
% 
% ; NAME:
% ;               bpass
% ; PURPOSE:
% ;               Implements a real-space bandpass filter that suppresses 
% ;               pixel noise and long-wavelength image variations while 
% ;               retaining information of a characteristic size.
% ;
% ; CATEGORY:
% ;               Image Processing
% ; CALLING SEQUENCE:
% ;               res = bpass( image_array, lnoise, lobject )
% ; INPUTS:
% ;               image:  The two-dimensional array to be filtered.
% ;               lnoise: Characteristic lengthscale of noise in pixels.
% ;                       Additive noise averaged over this length should
% ;                       vanish. May assume any positive floating value.
% ;               lobject: Integer length in pixels somewhat larger than
% ;                        a typical object.
% ; OUTPUTS:
% ;               res:    filtered image.
% ; PROCEDURE:
% ;               simple convolution yields spatial bandpass filtering.
% ; NOTES:
% ; MODIFICATION HISTORY:
% ;               Written by David G. Grier, The University of Chicago, 2/93.
% ;               Greatly revised version DGG 5/95.
% ;               Added /field keyword JCC 12/95.
% ;               Memory optimizations and fixed normalization, DGG 8/99.
%                 Converted to Matlab by D.Blair 4/2004-ish
%                 Fixed some bugs with conv2 to make sure the edges are
%                 removed D.B. 6/05
%                 Removed inadvertent image shift ERD 6/05
%                 Added threshold to output.  Now sets all pixels with
%                 negative values equal to zero.  Gets rid of ringing which
%                 was destroying sub-pixel accuracy, unless window size in
%                 cntrd was picked perfectly.  Now centrd gets sub-pixel
%                 accuracy much more robustly ERD 8/24/05

%                 Refactored for clarity and converted all convolutions to
%                 use column vector kernels for speed.  Running on my 
%                 macbook, the old version took ~1.3 seconds to do
%                 bpass(image_array,1,19) on a 1024 x 1024 image; this
%                 version takes roughly half that.
%                 JWM 6/07

% ;
% ;       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
% ;       David G. Grier.  It should be considered 'freeware'- and may be
% ;       distributed freely in its original form when properly attributed.
% 
%   

% Perform a bandpass by convolving with an appropriate kernel.  You can
% think of this as a two part process.  First, a lowpassed image is
% produced by convolving the original with a gaussian.  Next, a second
% lowpassed image is produced by convolving the original with a boxcar
% function. By subtracting the boxcar version from the gaussian version, we
% are using the boxcar version to perform a highpass.
% 
% original - lowpassed version of original => highpassed version of the
% original
% 
% Performing a lowpass and a highpass results in a bandpassed image.

  normalize = @(x) x/sum(x);

  image_array = double(image_array);
  r = -round(lobject):round(lobject);
  gaussian_kernel = normalize(exp(-(r/(2*lnoise)).^2))';
  boxcar_kernel = ones(length(r),1)/length(r);

% JWM: Do a 2D convolution with the kernels in two steps each.  It is possible
% to do the convolution in only one step per kernel with 
%
  % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
  % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
% 
% but for some reason, this is slow.  The whole operation could be reduced
% to a single step using the associative and distributive properties of
% convolution:
%
  % filtered = conv2(image_array,gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,'same');
%
% But this is also comparatively slow (though inexplicably faster than the
% above).  It turns out that convolving with a column vector is faster than
% convolving with a row vector, so instead of transposing the kernel, the
% image is transposed twice.

  gconv = conv2(image_array',gaussian_kernel,'same');
  gconv = conv2(gconv',gaussian_kernel,'same');

  bconv = conv2(image_array',boxcar_kernel,'same');
  bconv = conv2(bconv',boxcar_kernel,'same');

  filtered = gconv - bconv;

% JWM: The values on the borders of the result were convolved with 0s, so there
% will be fall off at the edges.  The old code set the affected pixels to
% zero, which seems unnecessary to me.  To reproduce the old results,
% uncomment the following 4 lines.
% ERD:  We zero them to signal that they're useless.

  filtered(1:(round(lobject)),:) = 0;
  filtered((end - lobject + 1):end,:) = 0;
  filtered(:,1:(round(lobject))) = 0;
  filtered(:,(end - lobject + 1):end) = 0;

% JWM: I question the value of zeroing out negative pixels.  It's a nonlinear
% operation which could potentially mess up our expectations about
% statistics.  Is there data on 'Now centroid gets subpixel accuracy much
% more robustly'?  To choose which approach to take, uncomment one of the
% following two lines.
% ERD: The negative values shift the peak if the center of the cntrd mask
% is not centered on the particle.

  % res = filtered;
  res = max(filtered,0);
