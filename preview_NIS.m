function hFigs = preview_NIS(filepath)
% Creates figures with having image from the first time/location point for 
% each channel in an NIS elements image file

if ~exist('filepath','var')
    [file,pth] = uigetfile(fullfile(last_dir,'*.nd2'),'Select NIS Elements image series file');
    if file==0
        return;
    end
    filepath = fullfile(file,pth);
end

if ~exist(filepath,'file')
    error('The specified file: %s does not exist',filepath);
end


%% Get file info
bfreader = bfGetReader(filepath);

numSeries = bfreader.getSeriesCount();
numChan = bfreader.getSizeC();

meta = bfreader.getMetadataStore();
bfreader.setSeries(0);
WIDTH = bfreader.getSizeX();
HEIGHT = bfreader.getSizeY();
numT = bfreader.getSizeT();


%% Load and display images
hFigs = gobjects(numChan,1);
for n=1:numChan
    hFigs(n) = figure;
    
    idx = bfreader.getIndex(0,n-1,0)+1;
    img = bfGetPlane(bfreader,idx);
    
    imagesc(img);
    axis image;
    colormap gray;
    
    title(sprintf('File: %s\nChannel: %d/%d  TimeStep: %d/%d Series: %d/%d',filepath,n,numChan,1,numT,1,numSeries),'Interpreter','none');
end

bfreader.close();
    