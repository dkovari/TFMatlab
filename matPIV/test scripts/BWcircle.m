function BW = BWcircle(r,loc,sz)
% Generates a size=[height,width] black & white image with of a circle of 
% radius r located at loc=[center_x, center_y].
% INPUT:
%   r: radius of circle
%   loc: [center_x, center_y] the location of the circle
%   size: [height, width] the size of the output image

[xx,yy] = meshgrid(1:sz(2),1:sz(1));

BW = sqrt( (xx-loc(1)).^2 + (yy-loc(2)).^2 )<= r;