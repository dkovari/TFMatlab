function clim = stackclim(imstack,method)

if nargin<2
    method = 'global';
end
clim=[0,0];
nFrames=size(imstack,3);
switch method
    case 'global'
        clim(1) = min(imstack(:));
        clim(2) = max(imstack(:));
    case 'average'
        Maxlist = nan(nFrames,1);
        Minlist = Maxlist;
        for f = 1:nFrames
            Maxlist(f) = nanmax(reshape(imstack(:,:,f),[],1));
            Minlist(f) = nanmin(reshape(imstack(:,:,f),[],1));
        end
        clim = [nanmean(Minlist),nanmean(Maxlist)];
    otherwise
        error('unknown method');
end