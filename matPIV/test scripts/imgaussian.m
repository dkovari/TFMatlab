function im = imgaussian(sig2,loc,sz)
% Generates a size=[height,width] grayscale image with of a gaussian with 
% varriance sig2 =[sig2x,sig2y] centered at loc=[center_x, center_y]. Image is normalized.
% INPUT:
%   sig2: [sig2x,sig2y] varriance
%   loc: [center_x, center_y] the location of the gaussian center
%   size: [height, width] the size of the output image

normalize2 = @(x) x/sum(sum(x));
[xx,yy] = meshgrid(1:sz(2),1:sz(1));

if numel(sig2)<2
    sig2=[sig2,sig2];
end


im = normalize2(exp(-((xx-loc(1) ).^2./(2*sig2(1)) + (yy-loc(2)).^2./(2*sig2(2)) )));