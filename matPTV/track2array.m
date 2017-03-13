function tarray = track2array(tracks)
% Convert output of track() to a more user-friendly t x [x,y] x ID array
%
% Input:
%   tracks - Output from tracks(), array of [x,y,t,id]
%       For now, t must be an integer corresponding to the frame number
%       t can start at either 0 or 1. If min(t)=1, it assumes that t=1 was
%       the first frame
%
% Output:
%   tarray: 3D array of track data
%       tarray(frame,xy,id) = [[x1,y1; x2,y2;...],...]

%Check to make sure time is an integer
if any(mod(tracks(:,3),1))
    error('time must be an integer, either starting at 0 or 1');
end

if min(tracks(:,3))==0
    tracks(:,3)=tracks(:,3)+1;
end
    
nF = max(tracks(:,3));
nID = max(tracks(:,4));

%Allocate empty array
tarray = NaN(nF,2,nID);

for r=1:size(tracks,1)
    tarray(tracks(r,3),:,tracks(r,4)) = tracks(r,1:2);
end
    



