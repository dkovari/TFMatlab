function OUT=shiftmat(IN,shiftvec,padval)
% out=shiftmat(in,shiftvec,padval)
% shifts the matrix IN while padding the shifted version by PADVAL.
% Inputs:
%   IN  input matrix
%   shiftvect [rows,cols]   amount to shift in give dimension
%   padval  Value to pad matrix
%       'circular','circ' simply use circshift

if nargin<3
    padval = 0;
end

if any(mod(shiftvec,1))
    error('shiftvec must contain only integers');
end

if(strcmpi(padval,'circular')||strcmpi(padval,'circ'))
    OUT=circshift(IN,shiftvec);
else
    if padval==0
        OUT = zeros(size(IN),'like',IN);
    elseif isnan(padval)
        OUT = NaN(size(IN),'like',IN);
    else
        OUT=padval*ones(size(IN),'like',IN); % Fill the output array with value to be padded.
    end
    
    if numel(shiftvec)>ndims(IN)
        error('shiftvec cannot have more dims than IN');
    end
    %build index array
    outidx = cell(numel(shiftvec),1);
    inidx = cell(numel(shiftvec),1);
    for d=1:numel(shiftvec)
        if shiftvec(d)<0
            outidx{d}=1:size(IN,d)+shiftvec(d);
            inidx{d} =-shiftvec(d)+1:size(IN,d);
        else
            outidx{d}= 1+shiftvec(d):size(IN,d);
            inidx{d} = 1:size(IN,d)-shiftvec(d);
        end
    end

     OUT(outidx{:})=IN(inidx{:});
end
