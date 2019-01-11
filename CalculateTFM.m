function TFMdata = CalculateTFM(varargin)
% Calculate Traction Force Data from ND2 Image files
%   Strain Field is computed using PTV 
%
% Optional Arguments:
%   'FilePath',string: Specify the ND2 file that should be loaded.
% 	'SeriesNum',int
% 	'CellChan',int
% 	'BeadChan',int
% 	'X_CROP',[int x1, int x2]
%     'Y_CROP',[int y1, int y2]
% 
% 	'Reference','PATH TO FILE' or 'same'
% 	'RefSeries',int
% 	'RefChan',int
% 	'RefFrame',int
% 
% 	'bpass_lnoise',int
% 	'bpass_sz',int
% 	'pkfnd_sz',int
% 	,'pkfnd_th',double
% 
% 	'YoungE',double
% 	'PoissonV',double
% 
% 	'SavePath','PATH TO OUTPUT'
% 	'SaveResults',true/false
% 
% 	'ForceMoviePath','PATH TO MOVIE'
% 
% 	'RelativeDisplacementOnly',false/true
% 	'MaxDisplacement',int
%
% 	'SaveSE',true/false: save plot of Strain Energy
% 	'CellImageCLim','average'/'global',[low,high]
% 	'MapLim','global'/'average',[low,high]
% 	'PlotStrain',true/false
%
% Output:
%   Unless directed not to save the data, the results are saved to a file
%   specified by 'SavePath' or the user file selection prompt. The
%   resulting matfile contains a struct:
%   TFMdata
%     TFMdata.Iref = Reference bead image used
%     TFMdata.bpass_lnoise = size of noise used for filtering
%     TFMdata.bpass_sz = particle size used for filtering
%     TFMdata.cnt_sz = particle size used for identifying particles
%     TFMdata.pkfnd_sz = particle size used for peak detection
%     TFMdata.pkfnd_th = intensity min used for peak detection
%     TFMdata.rF = reference frame number
%     TFMdata.nF = number of frames
%     TFMdata.RefFile = file used for reference
%     TFMdata.ImageFile = file containing original image data
%     TFMdata.SeriesNum = series used
%     TFMdata.BeadChan = channel used for bead images
%     TFMdata.CellChan = channel used for cell images
%     TFMdata.Ostack = cell images in a 3-d stack
%     TFMdata.imstack = bead images in a 3d stack
%     TFMdata.Bstack = filtered bead images in a 3d stack
%     TFMdata.Time = timepoints for each frame
%     TFMdata.Y_CROP = crop used on original data [1st,last]
%     TFMdata.X_CROP = crop used on original data [1st,last]
%     TFMdata.cnt = cell array listig parcle centers (only first cell
%                   contains data)
%     TFMdata.Vxx = location of displacement vectors along x
%     TFMdata.Vyy = location of displacement vectors along y
%     TFMdata.Vqx = x-value of displacement vectors
%     TFMdata.Vqy = y-value of displecement
%     TFMdata.SF = StressField;
%     TFMdata.dx = width of displacement/stress field windows in m
%     TFMdata.dy = height of displacement/stress field windows in m
%     TFMdata.dH = height of windows in pixels
%     TFMdata.dW = width of windoes in pixels
%     TFMdata.E = Youngs modulus of gel
%     TFMdata.v = Poisson's ration of gel
%     TFMdata.PX_SCALE = pixel scale (um/px)
%     TFMdata.SMAG = magnitude of stress for each window
%     TFMdata.SED = StrainEnergyDensity;
%     TFMdata.StrainEnergy = StrainEnergy, scalar total energy exerted
   

%% Parse Inputs
p = inputParser;
p.CaseSensitive = false;
addParameter(p,'FilePath',[]);
addParameter(p,'SeriesNum',[]);
addParameter(p,'CellChan',[]);
addParameter(p,'BeadChan',[]);
addParameter(p,'X_CROP',[]);
addParameter(p,'Y_CROP',[]);

addParameter(p,'Reference',[]); %file or "same"
addParameter(p,'RefSeries',[]);
addParameter(p,'RefChan',[]);
addParameter(p,'RefFrame',[]);

addParameter(p,'bpass_lnoise',[]);
addParameter(p,'bpass_sz',[]);
addParameter(p,'pkfnd_sz',[]);
addParameter(p,'pkfnd_th',[]);

addParameter(p,'YoungE',[]);
addParameter(p,'PoissonV',[]);

addParameter(p,'SavePath',[]);
addParameter(p,'SaveResults',true);

addParameter(p,'ForceMoviePath',[]);

addParameter(p,'RelativeDisplacementOnly',false);
addParameter(p,'MaxDisplacement',20);

addParameter(p,'SaveSE',true);
addParameter(p,'CellImageCLim','average');
addParameter(p,'MapLim','global');
addParameter(p,'PlotStrain',true);

parse(p,varargin{:});

%% Setup search path
if ~exist('piv_rec.m','file')   %check for matPIV library
    addpath(fullfile(fileparts(mfilename('fullpath')),'matPIV'));
    if ~exist('piv_rec.m','file')
        error('Could not find matPIV fitting library.');
    end
end
if ~exist('cntrd.m','file')   %check for matPTV library
    addpath(fullfile(fileparts(mfilename('fullpath')),'matPTV'));
    if ~exist('cntrd.m','file')
        error('Could not find matPTV fitting library.');
    end
end
if ~exist('find_tracks.m','file')   %check for matPTV library
    addpath(fullfile(fileparts(mfilename('fullpath')),'newTrack'));
    if ~exist('find_tracks.m','file')
        error('Could not find newTrack functions.');
    end
end

%% Load Bioformats Importer
if ~exist('bfGetReader.m','file') %check for bioformats_importer library
    addpath(fullfile(fileparts(mfilename('fullpath')),'bfmatlab'));
    if ~exist('bfGetReader.m','file')
        error('Could not find bfmatlab library.');
    end
end

%% Prompt User for data
persistent last_dir;

if ~isempty(p.Results.FilePath)
    if ~exist(p.Results.FilePath,'file')
        error('Specified File: %s does not exits',p.Results.FilePath);
    end
    
    [Dir,File,ext] = fileparts(p.Results.FilePath);
    File = [File,ext];
    
else
    %select file
    [File,Dir] = uigetfile(fullfile(last_dir,'*.nd2'),'Select NIS Elements image series file');
    if File==0
        return
    end
    if ~isempty(Dir)
        last_dir = Dir;
    end
end

%% Load File Info
hDlg = msgbox({'Loading ND2 File Info','Please wait'},'Loading...');
bfreader = bfGetReader(fullfile(Dir,File));
try
close(hDlg);
catch
end

numSeries = bfreader.getSeriesCount();
numChan = bfreader.getSizeC();

persistent SeriesNum;

%% Select Series and Channels

if ~isempty(p.Results.SeriesNum)
    SeriesNum = p.Results.SeriesNum;
end

if isempty(SeriesNum) || SeriesNum>numSeries
    SeriesNum = 1;
end


persistent BeadChan;

if ~isempty(p.Results.BeadChan)
    BeadChan = p.Results.BeadChan;
end

if isempty(BeadChan)||BeadChan>numChan
    BeadChan = 1;
end

persistent CellChan;



if numChan<2
    CellChan = [];
else
    if ~isempty(p.Results.CellChan)
        CellChan = p.Results.CellChan;
    elseif BeadChan~=1
        CellChan = 1;
    else
        CellChan = 2;
    end
end

first_loop = isempty(p.Results.SeriesNum)&&numSeries>1 || isempty(p.Results.BeadChan)&&numChan>1;
while first_loop||isnan(SeriesNum)||isnan(BeadChan)
    first_loop=false;
    if isnan(SeriesNum)
        SeriesNum = 1;
    end
    if isnan(BeadChan)
        BeadChan = 1;
    end
    
    if numSeries>1&&numChan>1
        prom = {sprintf('Position Series (1-%d)',numSeries);...
                sprintf('Bead Channel (1-%d)', numChan);...
                sprintf('Cell Channel (1-%d)', numChan)};
        defAns = {sprintf('%d',SeriesNum);sprintf('%d',BeadChan);sprintf('%d',CellChan)};
        answer = inputdlg(prom,'Select Series and Channel',1,defAns);
        
        if isempty(answer)
            bfreader.close()
            return;
        end
        
        SeriesNum = str2double(answer{1});
        BeadChan = str2double(answer{2});
        CellChan = str2double(answer{3});

    elseif numSeries>1
        prom = sprintf('Position Series (1-%d)',numSeries);
        defAns = {sprintf('%d',SeriesNum)};
        
        answer = inputdlg({prom},'Select Series (Assuming Bead=Chan1)',1,defAns);
        
        if isempty(answer)
            bfreader.close()
            return;
        end
        
        SeriesNum = str2double(answer{1});
    elseif numChan>1
        prom = {sprintf('Bead Channel (1-%d)', numChan);...
                sprintf('Cell Channel (1-%d)', numChan)};
        defAns = {sprintf('%d',BeadChan);sprintf('%d',CellChan)};
        answer = inputdlg(prom,'Select Channel (Only one Position)',1,defAns);
        if isempty(answer)
            bfreader.close()
            return;
        end
        
        BeadChan = str2double(answer{1});
        CellChan = str2double(answer{2});
    end
    if isempty(SeriesNum)||SeriesNum<1||SeriesNum>numSeries
        SeriesNum = NaN;
    end
    if isempty(BeadChan)||BeadChan<1||BeadChan>numChan
        BeadChan = NaN;
    end
end

%% Load Images
meta = bfreader.getMetadataStore();
bfreader.setSeries(SeriesNum-1);
WIDTH = bfreader.getSizeX();
HEIGHT = bfreader.getSizeY();
numT = bfreader.getSizeT();
PxScale = double(meta.getPixelsPhysicalSizeY(0).value());

USE_CELL_IMAGE = ~isempty(CellChan);

nF = numT;


Time = NaN(nF,1);

imstack = zeros(HEIGHT,WIDTH,nF);
if USE_CELL_IMAGE
    Ostack = zeros(HEIGHT,WIDTH,nF);
end

hWait = waitbar(0,'Loading Images');
for f=1:nF
    idx = bfreader.getIndex(0,BeadChan-1,f-1)+1;
    imstack(:,:,f) = bfGetPlane(bfreader,idx);
    
    %get timestamp
    dT = meta.getPlaneDeltaT(SeriesNum-1,idx-1);
    Time(f) = dT.value(ome.units.UNITS.S).doubleValue();
    
    if USE_CELL_IMAGE
        idx = bfreader.getIndex(0,CellChan-1,f-1)+1;
        Ostack(:,:,f) = bfGetPlane(bfreader,idx);
    end
    
    waitbar(f/nF,hWait);
end
bfreader.close();
try
close(hWait);
catch
end

%% Prompt for crop
if ~isempty(p.Results.X_CROP)
    if ischar(p.Results.X_CROP)&&strcmpi(p.Results.X_CROP,'full')
        X_CROP = [1,WIDTH];
    elseif p.Results.X_CROP(1)>0 && p.Results.X_CROP(1)<=WIDTH...
        && p.Results.X_CROP(2)>p.Results.X_CROP(1) && p.Results.X_CROP(2)<=WIDTH

        X_CROP = p.Results.X_CROP;
    end
end
if ~isempty(p.Results.Y_CROP)
    if ischar(p.Results.Y_CROP)&&strcmpi(p.Results.Y_CROP,'full')
        Y_CROP = [1,HEIGHT];
    elseif  p.Results.Y_CROP(1)>0 && p.Results.Y_CROP(1)<=HEIGHT...
        && p.Results.Y_CROP(2)>p.Results.Y_CROP(1) && p.Results.Y_CROP(2)<=HEIGHT
    
        Y_CROP = p.Resuts.Y_CROP;
    end
else
    answer = questdlg('Do you want to crop the images?');
    if strcmpi(answer,'cancel')
        return
    end
    if strcmpi(answer,'yes')
        if USE_CELL_IMAGE
            [Y_CROP,X_CROP] = uicropstack(Ostack,'colormap','gray','clim','average');
        else
            [Y_CROP,X_CROP] = uicropstack(imstack,'colormap','jet','clim','auto');
        end
        drawnow();



    else
        Y_CROP = [1,HEIGHT];
        X_CROP = [1,WIDTH];
    end
end

%apply crop
imstack = imstack(Y_CROP(1):Y_CROP(2),X_CROP(1):X_CROP(2),:);  
if USE_CELL_IMAGE
    Ostack = Ostack(Y_CROP(1):Y_CROP(2),X_CROP(1):X_CROP(2),:);
end

%% Reference Frame

if ~isempty(p.Results.Reference)
    if strcmpi(p.Results.Reference,'same')
        button = 'Choose Frame';
    else
        if ~exist(p.Results.Reference,'file')
            error('Specified Reference: %s does not exist',p.Results.Reference);
        end
        button = 'Alternate Image';
    end
else
    button = questdlg('What frame should be used as a reference image?','Reference Frame','First Frame','Choose Frame', 'Alternate Image','First Frame');
end

switch button
    case 'First Frame'
        rF = 1;
        Iref = imstack(:,:,1);
        RefFile = fullfile(Dir,File);
    case 'Choose Frame'
        if isempty(p.Results.RefFrame)
            rF = NaN;
        else
            if ischar(p.Results.RefFrame)
                switch lower(p.Results.RefFrame)
                    case 'first'
                        rF = 1;
                    case 'last'
                        rF = bfreader.getSizeT();
                end
            else
                rF = p.Results.RefFrame;
            end
        end
        while isnan(rF)
            answer = inputdlg(sprintf('Reference Frame (1-%d)',nF),'Reference',1,{'1'});
            if isempty(answer)
                return;
            end
            rF = str2double(answer);
            if rF<1||mod(rF,1)
                rF = NaN;
            end
        end
        Iref = imstack(:,:,rF);
        RefFile = fullfile(Dir,File);
    case 'Alternate Image'
        if ~isempty(p.Results.Reference)
            [RefDir,RefFile,ext] = fileparts(p.Results.Reference);
            RefFile = [RefFile,ext];
        else
            %% Prompt for File
            [RefFile,RefDir] = uigetfile(fullfile(last_dir,'*.nd2'),'Select NIS Elements image series file');
            if RefFile==0
                return;
            end
        end
        
        %% Load File Info
        hDlg = msgbox({'Loading ND2 File Info','Please wait'},'Loading...');
        bfreader = bfGetReader(fullfile(RefDir,RefFile));
        try
        close(hDlg);
        catch
        end

        RnumSeries = bfreader.getSeriesCount();
        RnumChan = bfreader.getSizeC();
        
        if RnumSeries == 1
            RSeriesNum = 1;
        else
            if ~isempty(p.Results.RefSeries)
                RSeriesNum = p.Results.RefSeries;
            else
                RSeriesNum = SeriesNum;
            end
        end
        
        if RnumChan==1
            RBeadChan = 1;
        else
            if ~isempty(p.Results.RefChan)
                RBeadChan = p.Results.RefChan;
            else
                RBeadChan = BeadChan;
            end
        end
        
        first_loop = isempty(p.Results.RefSeries)&&RnumSeries>1 ||...
                        isempty(p.Results.RefChan)&&RnumChan>1;
        while first_loop||isnan(RSeriesNum)||isnan(RBeadChan)
            first_loop = false;
            
            if RnumSeries>1&&RnumChan>1
                prom = {sprintf('Position Series (1-%d)',RnumSeries);...
                        sprintf('Bead Channel (1-%d)', RnumChan)};
                defAns = {sprintf('%d',RSeriesNum);sprintf('%d',RBeadChan)};
                answer = inputdlg(prom,'Select Series and Channel',1,defAns);

                if isempty(answer)
                    bfreader.close()
                    return;
                end

                RSeriesNum = str2double(answer{1});
                RBeadChan = str2double(answer{2});

            elseif RnumSeries>1
                prom = sprintf('Position Series (1-%d)',RnumSeries);
                defAns = {sprintf('%d',RSeriesNum)};

                answer = inputdlg({prom},'Select Series (Assuming Bead=Chan1)',1,defAns);

                if isempty(answer)
                    bfreader.close()
                    return;
                end

                RSeriesNum = str2double(answer{1});
            elseif RnumChan>1
                prom = {sprintf('Bead Channel (1-%d)', RnumChan)};
                defAns = {sprintf('%d',RBeadChan)};
                answer = inputdlg(prom,'Select Channel (Only one Position)',1,defAns);
                if isempty(answer)
                    bfreader.close()
                    return;
                end

                RBeadChan = str2double(answer{1});
            end
            if isempty(RSeriesNum)||RSeriesNum<1||RSeriesNum>numSeries
                RSeriesNum = NaN;
            end
            if isempty(RBeadChan)||RBeadChan<1||RBeadChan>numChan
                RBeadChan = NaN;
            end 
        end
        RefFile = fullfile(RefDir,RefFile);
        %% Select Frame
        bfreader.setSeries(RSeriesNum-1);
        
        if isempty(p.Results.RefFrame)
            if bfreader.getSizeT()==1
                rF = 1;
            else
                rF = NaN;
            end
        else
            if ischar(p.Results.RefFrame)
                switch lower(p.Results.RefFrame)
                    case 'first'
                        rF = 1;
                    case 'last'
                        rF = bfreader.getSizeT();
                end
            else
                rF = p.Results.RefFrame;
            end
        end
        while isnan(rF)
            answer = inputdlg(sprintf('Frame Number (1-%d)',bfreader.getSizeT()),'Reference Frame',1,{'1'});
            if isempty(answer)
                bfreader.close();
                return;
            end
            rF = str2double(answer);
            if rF<1||mod(rF,1)
                rF = NaN;
            end
        end
        
        idx = bfreader.getIndex(0,RBeadChan-1,rF-1)+1;
        
        %% Set Image
        Iref = bfGetPlane(bfreader,idx);
        bfreader.close();
        
        if any(size(Iref)~=[HEIGHT,WIDTH])
            error('Reference Image is not the same dimensions as uncropped sequence images');
        end
        Iref = Iref(Y_CROP(1):Y_CROP(2),X_CROP(1):X_CROP(2));
end

%% Particle Tracking Parameters
persistent bpass_lnoise;
persistent bpass_sz;
persistent pkfnd_sz;

if ~isempty(p.Results.bpass_lnoise)
    bpass_lnoise = p.Results.bpass_lnoise;
end

if isempty(bpass_lnoise)
    bpass_lnoise = 1;%approximate size of noise, in pixles (usually 1)
end

if ~isempty(p.Results.bpass_sz)
    bpass_sz = p.Results.bpass_sz;
end

if isempty(bpass_sz)
    bpass_sz = 7; %approximate particle size, in pixels. objects larger and smaller will be filtered with a spatial bandpass
end

if ~isempty(p.Results.pkfnd_sz)
    pkfnd_sz = p.Results.pkfnd_sz;
end
if isempty(pkfnd_sz)
    pkfnd_sz = 7;%approx. particle size, in pixels. size of the window used to find particle center of mass
end

first_loop = isempty(p.Results.bpass_lnoise)||isempty(p.Results.bpass_sz)||isempty(p.Results.pkfnd_sz);
while first_loop||isnan(bpass_lnoise)||isnan(bpass_sz)||isnan(pkfnd_sz)
    first_loop = false;
    
    answer = inputdlg(...
        {'Band Pass Noize size (px)';...
         'Band Pass Particle Size (px)';...
         'Peak-Find Window (px)'},...
         'Tracking Parameters',1,...
         {num2str(bpass_lnoise);...
         num2str(bpass_sz);...
         num2str(pkfnd_sz)});
     if isempty(answer)
         return;
     end
     bpass_lnoise = str2double(answer{1});
     bpass_sz = str2double(answer{2});
     pkfnd_sz = str2double(answer{3});
end      
persistent pkfnd_th;    %minimum peak intensity after bandpass operation

if ~isempty(p.Results.pkfnd_th)
    pkfnd_th = p.Results.pkfnd_th;
end

if isempty(pkfnd_th)
    pkfnd_th = 100;
end
cnt_sz = pkfnd_sz;         %particle size, in pixels

%% Calc Bpass stack
hWait = waitbar(0,'Filtering Image Noise');
Bstack = zeros(size(imstack));
for f=1:nF
    Bstack(:,:,f) = bpass(imstack(:,:,f),bpass_lnoise,bpass_sz);
    waitbar(f/nF,hWait);
end
try
    close(hWait);
catch
end

%% PIV drift correction
%Use coarse PIV to compensate for drift
scl = 3;
Iref = bpass(Iref,bpass_lnoise,bpass_sz);%imstack(:,:,rF);
Irefs = imresize(Iref,1/scl);
%imstack2=imstack;
driftXY = zeros(1,2,nF);
hWait = waitbar(0,'Coarse Drift Correction');
for f=1:nF
    %use xcorr to find approximate shift
    %[dXY,~] = imxcorr(Irefs,imresize(imstack(:,:,f),1/scl));
    [dXY,~] = imxcorr(Irefs,imresize(Bstack(:,:,f),1/scl));
    driftXY(1,:,f) = scl*dXY(:);
    %shift matrix, padd with zeros
    imstack(:,:,f) = shiftmat(imstack(:,:,f),-round([driftXY(1,2,f),driftXY(1,1,f)]),NaN);
    Bstack(:,:,f) = shiftmat(Bstack(:,:,f),-round([driftXY(1,2,f),driftXY(1,1,f)]),NaN);
    if USE_CELL_IMAGE
        Ostack(:,:,f) = shiftmat(Ostack(:,:,f),-round([driftXY(1,2,f),driftXY(1,1,f)]),'circ');
    end
    waitbar(f/nF,hWait);
end
try
close(hWait);
catch
end

%% Prompt for bead brightness
if isempty(p.Results.pkfnd_th)
    answer = questdlg('Do you want select bead brightness using image?');
    if strcmpi(answer,'yes')
        pkfnd_th = stackfig_threshold(Bstack,pkfnd_th);
    end
    val=NaN;
    while isnan(val)
        answer = inputdlg({'Brightness Threshold'},'Brightness',1,{num2str(pkfnd_th)});
        if ~isempty(answer)
            val = str2double(answer{1});
            if ~isnan(val)
                pkfnd_th = val;
            end
        else
            break
        end
    end
end

%% find Centroids
hWait = waitbar(0,'Finding Centroids');
cnt = cell(nF+1,1);
% first cnt will be the centroids in the refernce frame
PKS = pkfnd(Iref,pkfnd_th, pkfnd_sz);
c = cntrd(Iref,PKS,cnt_sz);
c(isnan(c(:,1)),:) = [];
cnt{1} = c(:,1:2);
% remaining cnt{} are the image frames
for f = 1:nF
    PKS = pkfnd(Bstack(:,:,f),pkfnd_th, pkfnd_sz);
    

    c = cntrd(Bstack(:,:,f),PKS,cnt_sz);
    c(isnan(c(:,1)),:) = [];
    
%     figure(1);
%     clf;
%     hold on;
%     plot(PKS(:,1),PKS(:,2),'o','markersize',20);
%     plot(c(:,1),c(:,2),'+','markersize',10);
%     title(sprintf('frame %d',f));
% %     pause;
%     drawnow;
    
    cnt{f+1} = c(:,1:2);
    waitbar(f/nF,hWait);
end
try
close(hWait);
catch
end

%% Track Particles
hDlg = msgbox({'Finding Particle Tracks','Please wait, this will take a while.'},'Tracking...');
if p.Results.RelativeDisplacementOnly
    new_tracks = find_displacement(cnt,1,'MaxDisp',p.Results.MaxDisplacement);
else
    new_tracks = find_tracks(cnt,'MaxDisp',p.Results.MaxDisplacement,'Memory',10,'ExpandingSearch',false);
end
try
close(hDlg);
catch
end

%ensure that each track is present at the reference frame
for t=numel(new_tracks):-1:1
    if isnan(new_tracks{t}(1,1))
        new_tracks(t) = [];
    end
end

%convert to 3d array
nID = numel(new_tracks);

if nID<1
    error('no tracks')
end

tracks = NaN(nF+1,2,nID);
for t=1:nID
    tracks(:,:,t) = new_tracks{t};
end

%% Calc Displacements
disptracks = bsxfun(@minus,tracks(2:end,:,:),tracks(1,:,:));

%calc average drift
meanDrift = nanmean(disptracks,3);
%correct for drift
tracks(2:end,:,:) = bsxfun(@minus,tracks(2:end,:,:),meanDrift);
disptracks = bsxfun(@minus,disptracks,meanDrift);

%% Bin Size
[H,W,~] = size(imstack);
% calculate average distance between particles
% H*W = dW*dH*size(tracks,3) one particle per bin
dW = sqrt(2*H*W/size(tracks,3)); % average of 2 beads per bin
%dW = 20;%round(1.5*max(W/sqrt(size(tracks,3)),H/sqrt(size(tracks,3))));
dH = dW;

%% Filter Displacements
% %cut-off factor
% bad_disp = false(size(disptracks,1),size(disptracks,3));
% hWait = waitbar(0,'Filtering Displacements');
% fact = 10;
% hoodfact = 3;
% for f=1:size(disptracks,1)
%     d = sqrt( sum(disptracks(f,:,:).^2,2));
%     medD = nanmedian(d)
%     idx = find(d>fact*medD);
%     for j = reshape(idx,1,[])
%         j
%         d(j)
%         %find neighbors in the same bin
%         nind = find(sum(bsxfun(@minus,tracks(1,:,:),tracks(1,:,j)).^2,2) < dW^2)
%         nind(nind==j) = [];
%         nind
%         hoodmean = nanmean(d(nind))
%         if ~isnan(hoodmean)&&d(j)>hoodmean*hoodfact
%             bad_disp(f,j) = true;
%         end
%     end
%     waitbar(f/nF,hWait);
% end
% try
% close(hWait);
% catch
% end

%% Interpolate Displacement
% Regularize displacement calculations
hWait = waitbar(0,'Interpolating Displacements');
%grid spacing in pixels
%dW = 16;
%dH = 16;


%Grid Points
[Vxx,Vyy] = meshgrid(0.5:dW:W+0.5,0.5:dH:H+0.5);

Vqx = zeros([size(Vxx),nF]);
Vqy = zeros([size(Vxx),nF]);

for f=1:nF
    %get rid of nans
    xyd = [reshape(tracks(1,1,:),[],1), ...
        reshape(tracks(1,2,:),[],1),...
        reshape(disptracks(f,1,:),[],1),...
        reshape(disptracks(f,2,:),[],1)];
    
    xyd(isnan(xyd(:,1)),:) = [];
    xyd(isnan(xyd(:,3)),:) = [];
    
    Fx = scatteredInterpolant(xyd(:,1),xyd(:,2),xyd(:,3));
    Fy = scatteredInterpolant(xyd(:,1),xyd(:,2),xyd(:,4));
    
    Vqx(:,:,f) = Fx(Vxx,Vyy);
    Vqy(:,:,f) = Fy(Vxx,Vyy);
    
    %fix edges to zero
    Vqx(:,[1,end],f) = 0;
    Vqx([1,end],:,f) = 0;
    Vqy(:,[1,end],f) = 0;
    Vqy([1,end],:,f) = 0;
    waitbar(f/nF,hWait);
end
try
close(hWait);
catch
end

%% Median filter displacement vectors
hWait = waitbar(0,'Filtering & Smoothing Displacement Vectors');

for f=1:nF
    if f==rF
        continue;
    end
    %filter by residual, with plainar fit to eliminate errors with PTV
    %[Vqx(:,:,f),Vqy(:,:,f),Vqxn(:,:,f),Vqyn(:,:,f)]=residfilt(Vqx(:,:,f),Vqy(:,:,f),2,'radius',2,'error',0.1,'method','linear');
    [Vqx(:,:,f),Vqy(:,:,f)]=residfilt(Vqx(:,:,f),Vqy(:,:,f),2,'radius',2,'error',0.1,'method','linear');
    Vqx(:,:,f)=wiener2(Vqx(:,:,f),[5,5]);   %wiener filter to smooth out noise
    Vqy(:,:,f)=wiener2(Vqy(:,:,f),[5,5]);   %wiener filter to smooth out noise
    waitbar(f/nF,hWait);
end
try
close(hWait);
catch
end

%% Prompt for Young's Modulous and Poisson's Ratio
persistent YoungE;
persistent PoissonV;

if ~isempty(p.Results.YoungE)
    YoungE = p.Results.YoungE;
end
if isempty(YoungE)
    YoungE = 7.2e3;%young's modulous;
end
if ~isempty(p.Results.PoissonV)
    PoissonV = p.Results.PoissonV;
end
if isempty(PoissonV)
    PoissonV = 0.46; %poisson's ratio
end

first_loop = isempty(p.Results.YoungE)||isempty(p.Results.PoissonV);
while first_loop||isnan(YoungE)||isnan(PoissonV)
    first_loop = false;
    
    answer = inputdlg(...
        {'Young''s Modulous (Pa)';...
         'Poisson''s Ratio'},...
         'Gel Stiffness',1,...
         {num2str(YoungE);...
         num2str(PoissonV)});
     if isempty(answer)
         return;
     end
     YoungE = str2double(answer{1});
     PoissonV = str2double(answer{2});
end 

%% Calculate Stress-Field

% Calculate stress using FTTC
hWait = waitbar(0,'Calculating Traction Force using FTTC');

%Force parameters
PX_SCALE = PxScale*1e-6;%0.107498 * 1e-6; %um/px

dx = PX_SCALE*dW;
dy = PX_SCALE*dH;
[Ny,Nx,~] = size(Vqx);
StressField = zeros(Ny,Nx,2,nF);
for f=1:nF
    [StressField(:,:,1,f),StressField(:,:,2,f)] = disp2stressFTTC(PX_SCALE*Vqx(:,:,f),PX_SCALE*Vqy(:,:,f),dx,dy,YoungE,PoissonV);
    waitbar(f/nF,hWait);
end
try
close(hWait);
catch
end

%% Calculate Stress mag
Stress_mag = zeros(size(Vqx));
nF = size(Vqx,3);
for f=1:nF
    Stress_mag(:,:,f) = sqrt(real(StressField(:,:,1,f)).^2+real(StressField(:,:,2,f)).^2);
end

%% Calculate Strain-energy density
StrainEnergyDensity = zeros(size(Vqx));
StrainEnergy = zeros(nF,1);
for f=1:nF
    StrainEnergyDensity(:,:,f) = PX_SCALE*0.5*real(StressField(:,:,1,f).*Vqx(:,:,f)+StressField(:,:,2,f).*Vqy(:,:,f));
    StrainEnergy(f) = sum(sum(StrainEnergyDensity(:,:,f)))*dx*dy;
end

%% Store Data in struct
if ~USE_CELL_IMAGE
    Ostack = [];
end

TFMdata = struct();
TFMdata.Iref = Iref;
TFMdata.bpass_lnoise = bpass_lnoise;
TFMdata.bpass_sz = bpass_sz;
TFMdata.cnt_sz = cnt_sz;
TFMdata.pkfnd_sz = pkfnd_sz;
TFMdata.pkfnd_th = pkfnd_th;
TFMdata.rF = rF;
TFMdata.nF = nF;
TFMdata.RefFile = RefFile;
TFMdata.ImageFile = fullfile(Dir,File);
TFMdata.SeriesNum = SeriesNum;
TFMdata.BeadChan = BeadChan;
TFMdata.CellChan = CellChan;
TFMdata.Ostack = Ostack;
TFMdata.imstack = imstack;
TFMdata.Bstack = Bstack;
TFMdata.Time = Time;
TFMdata.Y_CROP = Y_CROP;
TFMdata.X_CROP = X_CROP;
TFMdata.cnt = cnt;
TFMdata.driftXY = driftXY;
TFMdata.disptracks = disptracks;
TFMdata.tracks = tracks;
TFMdata.meanDrift = meanDrift;
TFMdata.Vxx = Vxx;
TFMdata.Vyy = Vyy;
TFMdata.Vqx = Vqx;
TFMdata.Vqy = Vqy;
TFMdata.SF = StressField;
TFMdata.dx = dx;
TFMdata.dy = dy;
TFMdata.dH = dH;
TFMdata.dW = dW;
TFMdata.E = YoungE;
TFMdata.v = PoissonV;
TFMdata.PX_SCALE = PX_SCALE;
TFMdata.SMAG = Stress_mag;
TFMdata.SED = StrainEnergyDensity;
TFMdata.StrainEnergy = StrainEnergy;

%% Save Data
if p.Results.SaveResults
    [~,name,~] = fileparts(File);
    if ~isempty(p.Results.SavePath)
        [dat_path,dat_file,ext] = fileparts(p.Results.SavePath);
        dat_file = [dat_file,ext];
    else
        [dat_file,dat_path] = uiputfile(fullfile(Dir,[name,sprintf('_XY%02d_TFMdata.mat',SeriesNum)]),'Save TFM Calculations?');
    end
    if dat_file~=0
        hDlg = msgbox({'Saving Data','Please wait, this will take a while.'},'Saving...');
        savefast(fullfile(dat_path,dat_file),...
            '-struct',TFMdata);
        try
        close(hDlg);
        catch
        end
    end   
end

%% Plot SE
if p.Results.SaveResults && p.Results.SaveSE
    hSE = figure();
    plot((TFMdata.Time-TFMdata.Time(1))/60,StrainEnergy,'-r');
    title('Strain Energy');
    ylabel('Strain Energy []');
    xlabel('Time [min]');
    
    [~,f,~] = fileparts(dat_file);
    saveas(hSE,fullfile(dat_path,[f,'_StrainEnergy.fig']));
    close(hSE);
end
%% View TFM Data
if isempty(p.Results.ForceMoviePath)
    TFMForceViewer(TFMdata,'CellImageCLim',p.Results.CellImageCLim,'MapLim',p.Results.MapLim,'PlotStrain',p.Results.PlotStrain,'FigureSize',[1080,720]);
else
    TFMForceViewer(TFMdata,'MoviePath',p.Results.ForceMoviePath,'CellImageCLim',p.Results.CellImageCLim,'MapLim',p.Results.MapLim,'PlotStrain',p.Results.PlotStrain,'FigureSize',[1080,720]);
end
%% Clear data if no output
if nargout<1
    clear TFMdata;
end