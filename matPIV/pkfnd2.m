function [PKS,IND] = pkfnd2(im,sz,thresh)
% Finds local maxima in an image, to pixel accuracy.
% Requirements:
%   MATLAB Image Processing Toolbox: imdilate
%       or
%   pkfnd(...) by E. Dufresne, et al. (included with this library)
% Inputs:
%   im: the grayscale image to find peaks in
%   sz: size of the local search window, either a number or MATLAB STREL
%       object
%   thresh: (optional) threshold
% Output:
%   PKS = [X,Y]
% =========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

if nargin<3
    thresh=mean(im(:));
end

try
    if isnumeric(sz)
        imt = imdilate(im,strel('square',sz));
    else
        imt = imdilate(im,sz);
    end
    imt(imt<thresh)=NaN;
    
    %pk = imt==im;
    
    
    IND = reshape(find(imt==im),[],1);
    %disp('old YX');
    %[Y,X] = ind2sub(size(im),find(imt==im));
%     disp('new Y')
    Y = rem(IND-1,size(im,1))+1;%
%     disp('new X');
    X = (IND-Y)/size(im,1) + 1;

    PKS = [X,Y];
catch
    %disp('using old pkfnd')
    [PKS,IND] = pkfnd(im,thresh,sz);
end

% figure(99);clf;
% hax1=subplot(1,2,1);
% imagesc('Parent',hax1,'CData',im); axis(hax1,'image'); colormap(hax1,'jet');
% hold(hax1,'on');
% plot(hax1,X,Y,'sw');
% hax2 = subplot(1,2,2);
% imt(isnan(imt)) = -10;
% imagesc('Parent',hax2,'CData',imt); axis(hax2,'image'); colormap(hax2,'jet');hold(hax2,'on');
% plot(hax2,X,Y,'sw');