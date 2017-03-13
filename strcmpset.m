function b = strcmpset(C1,C2)
%looks for strings in a set of strings
% Input:
%   C1: char array or cell of char arrays
%   C2: char array or cell of char arrays
% Output:
%   B = true/false, if C1 is a cell then size(B) = size(C1) with each
%   element specifying if C1{n} was matched to any string in C2

if iscell(C1)
    b = false(size(C1));
else
    
    if ~ischar(C1)
        error('C1 should be either a char array or a cell of char arrays');
    end
    b = any(strcmp(C1,C2));
    return
end

for n=1:numel(C1)
    b(n) = any(strcmp(C1{n},C2));
end