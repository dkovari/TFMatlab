function [x,y] = locPIV(imsz,sW,sH,eW,eH)
%Returns the location for velocity matrix returned by:
%   [...,x,y,...]=piv_rec(...,'startW,sW,'endW',eW,'startH',sH,'endH',eH)
%  Default:
%       startW = 64, startH = 64
%       endW = 16, endH = 16
% =========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-13: Initial creation (DTK)
%==========================================================================

if nargin<3
    sW = 64;
    sH = 64;
end
if nargin<5
    eW = 16;
    eH = 16;
end

dW = sW;
dH = sH;

x = 0.5:dW/2:imsz(2)+0.5;
y = 0.5:dH/2:imsz(1)+0.5;

if dW>eW
    dW = max(eW,dW/2);
    Wdone = false;
else
    Wdone = true;
end
if dH>eH
    dH = max(eH,dH/2);
    Hdone = false;
else
    Hdone = true;
end
while ~Wdone||~Hdone
    %subdivide positions
    x = reshape(x,1,[]);
    y = reshape(y,1,[]);
    x = reshape([x;[diff(x)/2,0]+x],1,[]); x(end)=[];
    y = reshape([y;[diff(y)/2,0]+y],1,[]); y(end)=[];

    %check if we can fit an extra row or col at end
    dx = x(end)-x(end-1);
    if x(end)+dx <= imsz(2)+0.5
        x = [x,x(end)+dx];
    end
    dy = y(end)-y(end-1);
    if y(end)+dy <= imsz(1)+0.5
        y = [y,y(end)+dy];
    end
    
    %reduce window size
    if dW>eW
        dW = max(eW,dW/2);
    else
        Wdone = true;
    end
    if dH>eH
        dH = max(eH,dH/2);
    else
        Hdone = true;
    end
end