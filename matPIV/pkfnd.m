function [out,ind]=pkfnd(im,th,sz,varargin)
% Finds local maxima in an image to pixel level accuracy.   
% This provides a rough guess of particle centers to be used by cntrd.m.
% Inspired by the lmx subroutine of Grier and Crocker's feature.pro
% INPUTS:
%   im: image to process, particle should be bright spots on dark
%       background with little noise, e.g. fluorescent images of µ-spheres
%       Use bpass.m to filter the image and eliminate noise

%   th: the minimum brightness of a pixel that might be local maxima. 
%       (NOTE: Make it big and the code runs faster
%       but you might miss some particles.  Make it small and you'll get
%       everything and it'll be slow.)
%   sz: (Optional) the target particle size. If two peaks are within sz/2
%       of each other only the brightest peak is preserved. It also removes
%       peaks within sz/2 of the boundary.
% OUTPUT:
%   out: [X,Y] and Nx2 array of the local maximum locations
%   ind: Nx1 array specifying the MATLAB index locations of each of the
%        peaks.  Useful for saving time when passing locations or values to
%        subsiquent functions.
%==========================================================================
%CREATED: Eric R. Dufresne, Yale University, Feb 4 2005

%% Change Log:
%MODIFIED: ERD, 5/2005, got rid of ind2rc.m to reduce overhead on tip by
%  Dan Blair;  added sz keyword 
% ERD, 6/2005: modified to work with one and zero peaks, removed automatic
%  normalization of image
% ERD, 6/2005: due to popular demand, altered output to give x and y
%  instead of row and column
% ERD, 8/24/2005: pkfnd now exits politely if there's nothing above
%  threshold instead of crashing rudely
% ERD, 6/14/2006: now exits politely if no maxima found
% ERD, 10/5/2006:  fixed bug that threw away particles with maxima
%  consisting of more than two adjacent points
% Dan Kovari, Curtis Lab Gatech, 02/2014: Clean up matlab syntax
% Dan Kovari (GT) 11/20/2014: Re-codded algorithm to use MATLAB
%   vectorization. Added the ind output variable.
%==========================================================================


%find all the pixels above threshold
%im=im./max(max(im)); 

p=inputParser;
p.CaseSensitive = false;
addParameter(p,'WarnBelowThresh',false);
parse(p,varargin{:});

if nargin<2
    th = mean(im(:));
end

ind=find(im > th);
[nr,nc]=size(im);
%tst=zeros(nr,nc);
if isempty(ind)
    out=[];
    if p.Results.WarnBelowThresh
        disp('nothing above threshold');
    end
    return;
end
%convert index from find to row and column
%[R,C] = ind2sub(size(im),ind);
R = rem(ind-1,size(im,1))+1;
C = (ind-R)/size(im,1) + 1;

%get rid of edge pixels
if nargin >=3   %user specified size, get rid of pixels too close to edge
    fd =find(R<floor(sz/2)+1|R>floor(nr-sz/2));
    ind(fd) = [];
    C(fd) = [];
    R(fd)=[];
    fd =find(C<floor(sz/2)+1|C>floor(nc-sz/2));
    ind(fd) = [];
    C(fd) = [];
    R(fd)=[];    
else
    fd =find(R<2|R>nr-1);
    ind(fd) = [];
    C(fd) = [];
    R(fd)=[];
    fd =find(C<2|C>nc-1);
    ind(fd) = [];
    C(fd) = [];
    R(fd)=[];  
end

% figure(99);
% imagesc(im); axis image; colormap jet;
% hold on;
% plot(C,R,'sw');

%iRC = ind;%R+(C-1)*nr;
iRCn = R+(C-2)*nr;
iRCp = R+(C)*nr;

iRnC = R-1+(C-1)*nr;
iRnCn = R-1+(C-2)*nr;
iRnCp = R-1+(C)*nr;

iRpC = R+1+(C-1)*nr;
iRpCn = R+1+(C-2)*nr;
iRpCp = R+1+(C)*nr;

ind = ind(...
            im(ind)>=im(iRCn)&...
            im(ind)>=im(iRCp)&...
            im(ind)>=im(iRnC)&...
            im(ind)>=im(iRnCn)&...
            im(ind)>=im(iRnCp)&...
            im(ind)>=im(iRpC)&...
            im(ind)>=im(iRpCn)&...
            im(ind)>=im(iRpCp)...
         );
     
ind = reshape(ind,[],1);

if nargin>=3&&~isempty(ind)    %user specified a min size range
    R = rem(ind-1,size(im,1))+1;
    C = (ind-R)/size(im,1) + 1;
    impk = zeros(size(im));
    impk(ind) = 1;
    
    for j=1:numel(ind)
        r1 = floor(R(j)-sz/2)+1;
        r2 = floor(R(j)+sz/2);
        c1 = floor(C(j)-sz/2)+1;
        c2 = floor(C(j)+sz/2);
        roi = impk( r1:r2,c1:c2);
        if sum(sum(roi))>1  %found a peak that was too close, choose birightest all other delete
            [~,mxi]=max(roi(:));
            roi = zeros(size(roi));
            roi(mxi)=1;
            impk(r1:r2,c1:c2) = roi;
        end
    end
    ind = reshape(find(impk),[],1);
end
if ~isempty(ind)
    R = rem(ind-1,nr)+1;
    C = (ind-R)/nr + 1;
    out = [C,R];
else
    out=[];
end

