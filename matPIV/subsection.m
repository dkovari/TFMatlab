function B = subsection(A,win,method)
% Subsection a matrix
% Input:
%   A input matrix
%   win [x1,y1,x2,y2] the subsection window corners
%   method: VAL (default = 0) value to pad edges
%           'mean' use mean value in subsection to pad edges
% Output:
%   B the sub-matrix specified by win
% =========================================================================
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Change Log:
%   2014-11-12: Initial creation (DTK)
%==========================================================================


if nargin<3
    method = 0;
end

sx = max(win(1),1);
sy = max(win(2),1);
ex = min(win(3),size(A,2));
ey = min(win(4),size(A,1));

B = A(sy:ey,sx:ex);

if all([sx,sy,ex,ey]==win)
    return;
end

if isnumeric(method)
    if any([sy-win(2),sx-win(1)])
        B = padarray(B,[sy-win(2),sx-win(1)],method,'pre');
    end
    if any ([win(4)-ey,win(3)-ex])
        B = padarray(B,[win(4)-ey,win(3)-ex],method,'post');
    end
elseif ischar(method)
    switch lower(method)
        case 'mean'
            method = mean(B(:));
            if any([sy-win(2),sx-win(1)])
                B = padarray(B,[sy-win(2),sx-win(1)],method,'pre');
            end
            if any ([win(4)-ey,win(3)-ex])
                B = padarray(B,[win(4)-ey,win(3)-ex],method,'post');
            end
        otherwise
            error('unknown method in in subsection');
    end
end

            