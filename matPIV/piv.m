function [vx2,vy2,x2,y2,vx_nan,vy_nan] = piv(im1,im2,dW,dH,varargin)
% PIV image tracking function
% Uses image cross-correlation to determine "flow-field" difference between
% two images.
% Inputs:
%   im1 first image
%   im2 second image
%   dW width (in pixels) of the cross-correlation window
%   dH height (in pixels) of the cross-correlation window
%       Without specifying any additional arguments, the flow field will be
%       calculated with windows which overlap by 50%.
% Optional inputs:
%   piv(...,x,y)    row vectors specifying the locations to use
%                   when compute correlation. 0.5 increments are
%                   interpreted as the space between two pixels.
%   piv(...,x,y,vx,vy)  specify initial guess for flow field
%               x,y are row vectors specifying locations of flow points
%               vx,vy are matricies of the flow vectors in x and y
%               directions, respectively.  The vectors should coorespond to
%               the locations given by [XX,YY] = meshgrid(x,y)
%   piv(...,x,y,vx,vy,'interp')
%               Use the specified flow field as above, but interpolate a
%               new mesh that is twice the density of the original.
%               Example:
%                   x=[1,2,3] => [1,1.5,2,2.5,3]
%   piv(...,x,y,vx,vy,'interp',method)
%               method: interpolation method to use [see interp2()]
%                   'linear' (default)
%                   'nearest'
%                   'cubic'
%                   'spline'
% Outputs:
%   vx2: array of flow vectors in x direction
%   vy2: array of flow vectors in the y direction
%   x2: vector listing x-location of vx2
%   y2: vector listing y-location of vy2
%   vx_nan: array of x-vectors, with NaNs indicating locations which
%           could not be computed because the cross-correlation did not
%           yield a peak of stron enough signal
%   vy_nan: array of y-vectors with NaNs...
% *************************************************************************
% Requirements:
%   findpeaks(...) included with this library
% =========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-13: Initial creation (DTK)
%==========================================================================
%tic
switch nargin
    case 4  %autocompute the locations of xcorr using window width/height
        x2 = 0.5:dW/2:size(im1,2)+0.5;
        y2 = 0.5:dH/2:size(im1,1)+0.5;
        vx2 = zeros(numel(y2),numel(x2));
        vy2 = zeros(numel(y2),numel(x2));
    case 6  %user specified xcorr locations
        x2= varargin{1}(:);
        y2 = varargin{2}(:);
        vx2 = zeros(numel(y2),numel(x2));
        vy2 = zeros(numel(y2),numel(x2));
    case 8  %user specified xcorr locations and initial flow
        x2= varargin{1}(:);
        y2 = varargin{2}(:);
        vx2 = varargin{3};
        vy2 = varargin{4};
        
        if any(size(vx2)~=[numel(y2),numel(x2)])||any(size(vy2)~=[numel(y2),numel(x2)])
            error('size of vx and vy must coorespond with [numel(y),numel(x)]');
        end
    case 9  %user specified that interpolation be used
        x2= varargin{1}(:);
        y2 = varargin{2}(:);
        vx2 = varargin{3};
        vy2 = varargin{4};
        
        if any(size(vx2)~=[numel(y2),numel(x2)])||any(size(vy2)~=[numel(y2),numel(x2)])
            error('size of vx and vy must coorespond with [numel(y),numel(x)]');
        end
        
        if (ischar(varargin{5})&&strcmpi(varargin{5},'interp')||...
                isscalar(varargin{5})&&varargin{5})
            x2 = reshape(x2,1,[]);
            y2 = reshape(y2,1,[]);
            x2 = reshape([x2;[diff(x2)/2,0]+x2],1,[]); x2(end)=[];
            y2 = reshape([y2;[diff(y2)/2,0]+y2],1,[]); y2(end)=[];
            
            vx2 = interp2(vx2,(1:0.5:size(vx2,2)),(1:0.5:size(vx2,1))');
            vy2 = interp2(vy2,(1:0.5:size(vy2,2)),(1:0.5:size(vy2,1))');
            
            %check if we can fit an extra row or col at end
            dx = x2(end)-x2(end-1);
            if x2(end)+dx <= size(im1,2)+0.5
                x2 = [x2,x2(end)+dx];
                vx2 = cat(2,vx2,vx2(:,end));    %pad end+1 with end
                vy2 = cat(2,vy2,vy2(:,end));    %pad end+1 with end
            end
            dy = y2(end)-y2(end-1);
            if y2(end)+dy <= size(im1,1)+0.5
                y2 = [y2,y2(end)+dy];
                vx2 = cat(1,vx2,vx2(end,:));    %pad end+1 with end
                vy2 = cat(1,vy2,vy2(end,:));    %pad end+1 with end
            end
        end
            
    case 10 %user specified method of intrpolation
        x2= varargin{1}(:);
        y2 = varargin{2}(:);
        vx2 = varargin{3};
        vy2 = varargin{4};
        
        if any(size(vx2)~=[numel(y2),numel(x2)])||any(size(vy2)~=[numel(y2),numel(x2)])
            error('size of vx and vy must coorespond with [numel(y),numel(x)]');
        end
        
        if (ischar(varargin{5})&&strcmpi(varargin{5},'interp')||...
                isscalar(varargin{5})&&varargin{5})
            x2 = reshape(x2,1,[]);
            y2 = reshape(y2,1,[]);
            x2 = reshape([x2;[diff(x2)/2,0]+x2],1,[]); x2(end)=[];
            y2 = reshape([y2;[diff(y2)/2,0]+y2],1,[]); y2(end)=[];
            
            vx2 = interp2(vx2,(1:0.5:size(vx2,2)),(1:0.5:size(vx2,1))',varargin{6});
            vy2 = interp2(vy2,(1:0.5:size(vy2,2)),(1:0.5:size(vy2,1))',varargin{6});
        end
        
        %check if we can fit an extra row or col at end
        dx = x2(end)-x2(end-1);
        if x2(end)+dx <= size(im1,2)+0.5
            x2 = [x2,x2(end)+dx];
            vx2 = cat(2,vx2,vx2(:,end));    %pad end+1 with end
            vy2 = cat(2,vy2,vy2(:,end));    %pad end+1 with end
        end
        dy = y2(end)-y2(end-1);
        if y2(end)+dy <= size(im1,1)+0.5
            y2 = [y2,y2(end)+dy];
            vx2 = cat(1,vx2,vx2(end,:));    %pad end+1 with end
            vy2 = cat(1,vy2,vy2(end,:));    %pad end+1 with end
        end
    otherwise
        error('incorrect number of arguments');
end

if nargout>4    %only create the nan array if user wants it
    vx_nan=NaN(size(vx2));
    vy_nan = NaN(size(vy2));
end

[xx2,yy2] = meshgrid(x2,y2);

%fprintf('running piv dW=%d dH=%d numXin=%d\n',dW,dH,numel(xx2));

for j=1:numel(vx2) %this could be parfor
    %calc new window limits
    win1 = [floor(xx2(j)-dW/2)+1,...
            floor(yy2(j)-dH/2)+1,...
            floor(xx2(j)+dW/2),...
            floor(yy2(j)+dH/2)];
    win2 = win1+round([vx2(j),vy2(j),vx2(j),vy2(j)]);
    
    %if either window is not completely in image, dont process
    if win1(1)>size(im1,2) || win1(3)>size(im1,2) || ...
            win1(2)>size(im1,1)||win1(4)>size(im1,1)||...
            win1(1)<1||win1(3)<1||...
            win1(2)<1||win1(4)<1||...
            win2(1)>size(im2,2)||win2(3)>size(im2,2)||...
            win2(2)>size(im2,1)||win2(4)>size(im2,1)||...
            win2(1)<1||win2(3)<1||...
            win2(2)<1||win2(4)<1
        dXY = [];
        SNR = 0;
    else  
        ext1 = win1;
        ext2 = win2;
        
        im1t = im1(ext1(2):ext1(4),ext1(1):ext1(3));
        im2t = im2(ext2(2):ext2(4),ext2(1):ext2(3));
        
        [dXY,SNR] = imxcorr(im1t,im2t);

    end
    
    if ~isempty(dXY)&&SNR>1.5
        vx2(j)=vx2(j)+dXY(1);
        vy2(j)=vy2(j)+dXY(2);
        if nargout>4    %only create the nan array if user wants it
            vx_nan(j) = vx2(j);
            vy_nan(j) = vy2(j);
        end
    end
end
