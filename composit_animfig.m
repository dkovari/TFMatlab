function [hFig,hAx,hImages,hCB] = composit_animfig(FrameData,varargin)
% Display a multi-frame "movie" of composit images
% The resulting gui resembles the GUI created by stackfig(...), which
% diplays the frames of a movie stored as a simple HxWx(RGBx)Frames movie.
% This view will accept an unlimited number of "movies" and overlay them.
% AlphaData for each frame is used to layer the images.
%
% Syntax::
%   composit_animfig(Mov1, Mov2, Mov3,...,Name,Value)
% Inputs
%   Mov1:
%       Images to be displayed as frames of the movie. Data can be numeic
%       arrays, cell arrays, or an animation structure
%    Numeric:
%       Mov1= H x W x nFrames array
%           each frame of the image will be interpreted as a HxW grayscale
%           image. The colormap and CLims are determined by the values of
%           the "colormap" and "CLim" parameters
%
%       Mov= H x W x 3 x nFrames
%           each frame is interpreted as an HxW RGB image
%
%       Mov1 = animation structure array (see image properties)
%          Mov1(f).CData = Image color data for frame f.
%                          If data is a matrix then colormap is used to
%                          convert it to RGB.
%              .CDataMapping = how image should be scaled with colormap  
%               .AlphaData = Transparency data for the image
%               .AlphaDataMapping = How transparency scales
%               .XData = [1x2] location image should be displayed along x
%               .YData = [1x2] location along y
%               .colormap = colormap that should be used for displaying
%                           image
% Parameters:
%   'colormap',... set colormap (see colormap() )
%   'Figure',hFig:  specify a figure to use
%   'ShowColorbar',true/false: show a customizable colorbar on right
%   'CustomColorbar',{Name,Value} specify colorbar axes properties
%   'ColorbarWidth',## width of colorbar in characters
%
%   'AlphaData',# default value for alphadata for images
%   'AlphaDataMapping','...': how alphadata should be scaled (see image)
% clim: default CLim to use for images that don't have CDataMapping Fields
%   'direct' = interpret as indexed array
%   'scaled' = rescale each frame to match colormap limits (like imagesc)
%   'global' = scale colormap to min/max of all images in stack
%   'average' = scale colormapt to average min/max of each image
%   [min,max] = specify specific value for clim
% user callback: 'frameupdate_fn',@(hFig,hAx,currentFrame)
%   Specify a callback function that should be executed when frame is
%   changed
%   access data using getappdata(hFig,...)
%  AppData Fields:
%       hImages
%       hAx
%       curFrame
%       nFrames
%       FrameData
%       nImages
%       hSld_Frame
%       hEdt_Frame
%       hTxt_Frame
%       frameupdate_fn
%       ExecuteMovie_Fn(filepath,varargin): function handle to save a movie
%           varargin parameters:
%               'Quality',default=75    the compression quality
%               'FrameRate',default=5   video frame rate
%   Additional variables can be included in the appdata using the 'appdata'
%   paramter
%
% Output:
%   hFig: handle to figure
%   hAx: handle to main axes
%   hImages: handle to images being shown on hAx
%   hCB: handle to customizeable colorbar (if ShowColorbar=true)
%% Copyright 2016 Daniel T. Kovari, Emory University
% All rights reserved.
%% Change Log:
%   2016-12-13: DTK
%       Initial File Creation

%% Input Parameters
p = inputParser();
p.CaseSensitive = false;
addParameter(p,'colormap','gray',@(x) ischar(x)||isnumeric(x)&&(size(x,2)==3));
addParameter(p,'CLim','direct',@(x) ischar(x)&&any(strcmpi(x,{'direct','scaled','global','average','manual'}))||isnumeric(x)&&numel(x)==2);
addParameter(p,'frameupdate_fn',[],@(x) isempty(x)||isa(x,'function_handle'));
addParameter(p,'Figure',[],@(x) isempty(x)||ishghandle(x));
addParameter(p,'appdata',struct(),@isstruct);
addParameter(p,'AlphaData',1);
addParameter(p,'AlphaDataMapping','none');
addParameter(p,'CustomColorbar',{});
addParameter(p,'ShowColorbar',false,@(x) islogical(x)&&isscalar(x));
addParameter(p,'ColorbarWidth',12,@(x) isscalar(x)&&isnumeric(x));
addParameter(p,'PreSaveFcn',[]);
addParameter(p,'PostSaveFcn',[]);

%% Validate/convert inputs
FDnames = {'CData','colormap','CDataMapping','XData','YData','AlphaData','AlphaDataMapping'};

FrameData = {FrameData};
% find last index of data
for k=1:numel(varargin)
   if ischar(varargin{k})
       k=k-1;
       break;
   end
end
FrameData = [FrameData,varargin{1:k}];

%% parse name,value pairs
parse(p,varargin{k+1:end});
defaultCmap = p.Results.colormap;

if ~isempty(p.Results.Figure)
    hFig = figure(p.Results.Figure);
    clf(hFig);
else
    hFig = figure();
end

defaultCDataMapping = p.Results.CLim;
% if isnumeric(p.Results.CLim)
%     defaultCDataMapping = 'scaled';
% else
%     defaultCDataMapping = p.Results.CLim;
% end

defaultAlphaData = p.Results.AlphaData;
defaultAlphaDataMapping = p.Results.AlphaDataMapping;

%% Assemble Frame Data
for fr = 1:numel(FrameData)
    %% n-d array data
    if isnumeric(FrameData{fr})
        FrameStruct = struct('CData',{});
        if ndims(FrameData{fr})==3
            for n=1:size(FrameData{fr},3)
                FrameStruct(n).CData = FrameData{fr}(:,:,n);
            end
        elseif ndims(FrameData{fr})==4 && size(FrameData{fr})==3
            for n=1:size(FrameData{fr},4)
                FrameStruct(n).CData = FrameData{fr}(:,:,:,n);
            end
        else
            error('Frame %d is numeric but is not HxWxN nor HxWx3xN',fr);
        end
        FrameData{fr} = FrameStruct;
    end
    %% Data is a cell array of images, convert to struct
    if iscell(FrameData{fr})
        %cdata
        [FrameStruct(1:numel(FrameData{fr})).CData] = deal(FrameData{fr}{:});
        FrameData{fr} = FrameStruct;
    end
    %% Data is already a struct
    if isstruct(FrameData{fr})
        %rename fields to match case
        fnames = fieldnames(FrameData{fr});
        for f=1:numel(fnames)
            ind = find(strcmpi(fnames{f},FDnames));
            if ~isempty(ind)
                if strcmp(fnames{f},FDnames{ind})
                    continue;
                end
                [FrameData{fr}.(FDnames{ind})] = FrameData{fr}.(fnames{f});
                FrameData{fr} = rmfield(FrameData{fr},fnames{f});
            end
        end

        if ~isfield(FrameData{fr},'CData')
            error('Input: frames struct must contain "CData" field');
        end
        %% add defaults
        if ~isfield(FrameData{fr},'CDataMapping')
            [FrameData{fr}.CDataMapping] = deal(defaultCDataMapping);
        end
        if ~isfield(FrameData{fr},'colormap')
            [FrameData{fr}.colormap] = deal(defaultCmap);
        end
        if ~isfield(FrameData{fr},'XData')
            for k=1:numel(FrameData{fr})
                FrameData{fr}(k).XData = [1,size(FrameData{fr}(k).CData,2)];
            end
        end
        if ~isfield(FrameData{fr},'YData')
            for k=1:numel(FrameData{fr})
                FrameData{fr}(k).YData = [1,size(FrameData{fr}(k).CData,1)];
            end
        end
        if ~isfield(FrameData{fr},'AlphaData')
            [FrameData{fr}.AlphaDataMapping] = deal(defaultAlphaDataMapping);
            for k=1:numel(FrameData{fr})
                FrameData{fr}(k).AlphaData = ~isnan(FrameData{fr}(k).CData(:,:,1));
                if all(FrameData{fr}(k).AlphaData)
                    FrameData{fr}(k).AlphaData = defaultAlphaData;
                end
            end
        end
        if ~isfield(FrameData{fr},'AlphaDataMapping')
            [FrameData{fr}.AlphaDataMapping] = deal(defaultAlphaDataMapping);
        end
    end
end

nFrames = 0;
for fr = 1:numel(FrameData)
    nFrames = max(nFrames,numel(FrameData{fr}));
end



%% Layout Figure
set(hFig,'units','characters');
fpos = get(hFig,'position');

%% Custom Colorbar
ColorbarParams = p.Results.CustomColorbar;
cbWidth = p.Results.ColorbarWidth;
UseCB = p.Results.ShowColorbar;
hCB = [];
if UseCB||~isempty(ColorbarParams)
    UseCB = true;
    hCB = axes('parent',hFig,...
        'units','characters',...
        'outerposition',[fpos(3)-cbWidth-2,5,cbWidth,fpos(4)-7],...
        'XTickLabel','',...
        'XTickLabelMode','manual',...
        'XTick',[],...
        'YAxisLocation','right',...
        'XColor','none',...
        ColorbarParams{:});
    hAxPos = [0,2,fpos(3)-cbWidth-2,fpos(4)-2];
else
    hAxPos = [0,2,fpos(3),fpos(4)-2];
end

hAx = axes('parent',hFig,...
            'units','characters',...
            'outerposition',hAxPos);

hSld_Frame = uicontrol('Parent', hFig,'Style','slider',...
         'Units','characters',...
         'Position',[1, 0, fpos(3)-15-1, 2],...
         'Min',1,...
         'Max',nFrames,...
         'Value',1,...
         'SliderStep',[1/nFrames,10/nFrames],...
         'Tag','hSld_Frame',...
         'Callback',@SldFrame_callback);

hEdt_Frame = uicontrol('Parent',hFig,...
        'Style','edit',...
        'units','characters',...
        'Position',[fpos(3)-15+1,0.7,7,1.6],...
        'String','1',...
        'HorizontalAlignment','center',...
        'Callback',@EditFrame_callback);
hTxt_Frame = uicontrol('Parent',hFig,...
    'Style','text',...
    'units','characters',...
    'position',[fpos(3)-15+1+7+.5,0,5,2],...
    'HorizontalAlignment','left',...
    'string',sprintf('/%d',nFrames));

%% Precalc ImageIntensities incase we need to use colormap
for f = 1:numel(FrameData)
    for n=1:numel(FrameData{f})
        FrameData{f}(n).Ilow = NaN;
        FrameData{f}(n).Ihigh = NaN;
        if ndims(FrameData{f}(n).CData)<3
            FrameData{f}(n).Ilow = nanmin(FrameData{f}(n).CData(:));
            FrameData{f}(n).Ihigh = nanmax(FrameData{f}(n).CData(:));
        end
    end
end
%% Convert colormapped images to rgb
for f = 1:numel(FrameData)
    for n=1:numel(FrameData{f})
        if ndims(FrameData{f}(n).CData)<3
            if ischar(FrameData{f}(n).colormap)
                cmap = eval([FrameData{f}(n).colormap,'(65536)']);
            else
                cmap = FrameData{f}(n).colormap;
            end
            if ischar(FrameData{f}(n).CDataMapping)
                switch lower(FrameData{f}(n).CDataMapping)
                    case 'direct'
                        %do nothing
                    case 'manual'
                        %do nothing
                    case 'scaled'
                        FrameData{f}(n).CDataMapping = [FrameData{f}(n).Ilow,FrameData{f}(n).Ihigh];
                    case 'global'
                        FrameData{f}(n).CDataMapping = [nanmin([FrameData{f}.Ilow]),nanmax([FrameData{f}.Ihigh])];
                    case 'average'
                        FrameData{f}(n).CDataMapping = [nanmean([FrameData{f}.Ilow]),nanmean([FrameData{f}.Ihigh])];
                end
            end
            %convert to RGB for all except first image
            if f~=1
                if ischar(FrameData{f}(n).CDataMapping) %still direct
                    FrameData{f}(n).CData = ind2rgb(FrameData{f}(n).CData,cmap);
                else %use numeric clim
                    FrameData{f}(n).CData = ind2rgb(gray2ind(mat2gray(FrameData{f}(n).CData,FrameData{f}(n).CDataMapping),size(cmap,1)),cmap);
                end
            end
        end
    end
    if f~=1
        FrameData{f} = rmfield(FrameData{f},{'colormap','CDataMapping','Ilow','Ihigh'});
    end
end


%% Draw Images
nImages = numel(FrameData);
hImages = gobjects(nImages,1);
hold(hAx,'on');

%first image is special case since it is not necessarily RGB
hImages(1) = image(hAx,'CData',FrameData{1}(1).CData,...
    'XData',FrameData{1}(1).XData,...
    'YData',FrameData{1}(1).YData,...
    'AlphaData',FrameData{1}(1).AlphaData,...
    'AlphaDataMapping',FrameData{1}(1).AlphaDataMapping,...
    'CDataMapping','scaled');
if ~ischar(FrameData{1}(1).CDataMapping)
    set(hAx,'CLim',FrameData{1}(1).CDataMapping);
end
colormap(hAx,FrameData{1}(1).colormap);

for n=2:nImages
    hImages(n) = image(hAx,'CData',FrameData{n}(1).CData);
    set(hImages(n),FrameData{n}(1));
end
axis(hAx,'image');
set(hAx,'ydir','reverse');

%% Set appdata for fig
appdata = p.Results.appdata;
appdata.hImages = hImages;
appdata.hAx = hAx;
appdata.curFrame = 1;
appdata.nFrames = nFrames;
appdata.FrameData = FrameData;
appdata.nImages = nImages;
appdata.hSld_Frame = hSld_Frame;
appdata.hEdt_Frame = hEdt_Frame;
appdata.hTxt_Frame = hTxt_Frame;
appdata.frameupdate_fn = p.Results.frameupdate_fn;
appdata.hCB = hCB;
appdata.cbWidth = cbWidth;

fld = fieldnames(appdata);
for n=1:numel(fld)
    setappdata(hFig,fld{n},appdata.(fld{n}));
end

set(hFig,'currentaxes',hAx);


set(hFig,'SizeChangedFcn',@ReSzFig);

%% Add Animation Menu
uimenu(hFig,'Label','Save Movie','Callback',@(~,~) AnimateFigure(hFig));

setappdata(hFig,'PreSaveFcn',p.Results.PreSaveFcn);
setappdata(hFig,'PostSaveFcn',p.Results.PostSaveFcn);

ReSzFig(hFig);

setappdata(hFig,'ExecuteMovie_Fn',@(fp,varargin) AnimateFigure(hFig,fp,varargin));
end

%% Resize Figure Function
function ReSzFig(hFig,~)
    orgUnits = get(hFig,'units');
    set(hFig,'units','characters');
    fpos = get(hFig,'position');
    hAx = getappdata(hFig,'hAx');
    hCB = getappdata(hFig,'hCB');
    cbWidth = getappdata(hFig,'cbWidth');
    
    hSld = getappdata(hFig,'hSld_Frame');
    hEdt = getappdata(hFig,'hEdt_Frame');
    hTxt = getappdata(hFig,'hTxt_Frame');
    
    HasCtrls = ishandle(hSld)&&ishandle(hEdt)&&ishandle(hTxt);
    
    
    if HasCtrls
        if ishghandle(hCB)
            set(hAx,'units','characters','outerposition',[0,2,fpos(3),fpos(4)-2]);
            %oAx = hAx.OuterPosition
            pAx = hAx.Position;
            tiAx = hAx.TightInset;
            
            if fpos(3)-(pAx(1)+pAx(3)+tiAx(3)) < cbWidth
                %'here'
                set(hCB,'units','characters','outerposition',[fpos(3)-cbWidth-2,5,cbWidth,fpos(4)-8]);
                set(hAx,'units','characters','outerposition',[0,2,fpos(3)-cbWidth-2,fpos(4)-2]);
            else
                set(hCB,'units','characters','outerposition',[pAx(1)+pAx(3)+tiAx(3),5,cbWidth,fpos(4)-8]);
            end
        else
            set(hAx,'units','characters','outerposition',[0,2,fpos(3),fpos(4)-2]);
        end
        try
        set(getappdata(hFig,'hSld_Frame'),'units','characters','Position',[1, 0, fpos(3)-15-1, 2]);
        set(getappdata(hFig,'hEdt_Frame'),'units','characters','Position',[fpos(3)-15+1,0.7,7,1.6]);
        set(getappdata(hFig,'hTxt_Frame'),'units','characters','position',[fpos(3)-15+1+7+.5,0,5,2]);
        catch
        end
    else
        if ishghandle(hCB)
            set(hAx,'units','characters','outerposition',[0,0,fpos(3),fpos(4)]);
            pAx = hAx.Position;
            tiAx = hAx.TightInset;
            
            if fpos(3)-(pAx(1)+pAx(3)+tiAx(3)) < cbWidth
                set(hCB,'units','characters','outerposition',[fpos(3)-cbWidth-2,3,cbWidth,fpos(4)-6]);
                set(hAx,'units','characters','outerposition',[0,0,fpos(3)-cbWidth-2,fpos(4)]);
            else
                set(hCB,'units','characters','outerposition',[pAx(1)+pAx(3)+tiAx(3),4,cbWidth,fpos(4)-6]);
            end
        else
            set(hAx,'units','characters','outerposition',[0,0,fpos(3),fpos(4)]);
        end
    end
    set(hFig,'units',orgUnits);
end
function SldFrame_callback(hObject,~)
hFig = get(hObject,'parent');
cF = round(get(hObject,'value'));
change_frame(hFig,cF);
end
    
function EditFrame_callback(hObject,~)
hFig = get(hObject,'parent');
cF = round(str2double(get(hObject,'string')));
if isnan(cF) || cF<1 || cF>getappdata(hFig,'nFrames')
    cF = getappdata(hFig,'curFrame');
    set(hObject,'string',num2str(cF));
    return;
end
change_frame(hFig,cF);

end

function change_frame(hFig,cF)
hAx = getappdata(hFig,'hAx');
%% Draw Images
nImages = getappdata(hFig,'nImages');
hImages = getappdata(hFig,'hImages');
FrameData = getappdata(hFig,'FrameData');
hold(hAx,'on');
for n=1:nImages
    if numel(FrameData{n})>=cF
        if n==1
            set(hImages(n),'CData',FrameData{n}(cF).CData,...
                'XData',FrameData{n}(cF).XData,...
                'YData',FrameData{n}(cF).YData,...
                'AlphaData',FrameData{n}(cF).AlphaData,...
                'AlphaDataMapping',FrameData{n}(cF).AlphaDataMapping,...
                'CDataMapping','scaled');
            colormap(hAx,FrameData{n}(cF).colormap);
            if ~ischar(FrameData{n}(cF).CDataMapping)
                %'here'
                %disp(FrameData{n}(cF).CDataMapping)
                %figure(2);imagesc(FrameData{n}(cF).CData)
                set(hAx,'CLim',FrameData{n}(cF).CDataMapping);
                %disp(hAx.CLim)
            end
        else
            set(hImages(n),FrameData{n}(cF));
        end
        set(hImages(n),'Visible','on');
    else
        set(hImages(n),'Visible','off');
    end
end
setappdata(hFig,'curFrame',cF);

%% Set slider and edit
hSld = getappdata(hFig,'hSld_Frame');
if ishandle(hSld)
    set(hSld,'value',cF);
end
hEdt = getappdata(hFig,'hEdt_Frame');
if ishandle(hEdt)
    set(hEdt,'string',num2str(cF));
end
%run user callback
fcn = getappdata(hFig,'frameupdate_fn');
if ~isempty(fcn)
    fcn(hFig,hAx,cF);
end
end

function AnimateFigure(hFig, movfilepath,varargin)

p=inputParser();
p.CaseSensitive = false;
addParameter(p,'Quality',75,@(x) isnumeric(x)&&isscalar(x)&&x>=0&&x<=100);
addParameter(p,'FrameRate',5,@(x) isnumeric(x)&&isscalar(x));
parse(p);

%% disable resize
fRsz = get(hFig,'Resize');
set(hFig,'Resize','off');

%% get appdata, apply presave fcn
PreSaveFcn = getappdata(hFig,'PreSaveFcn');

if ~isempty(PreSaveFcn)
    PreSaveFcn();
end

if nargin<2
    [mov_file, mov_path] = uiputfile('*.mp4','Save Animation?');

    if mov_file==0
        return;
    end
else
    [mov_path,mov_file,ext] = fileparts(movfilepath);
    mov_file = [mov_file,ext];
end

%% Delete Slider and frame counter
delete(getappdata(hFig,'hSld_Frame'));
delete(getappdata(hFig,'hEdt_Frame'));
delete(getappdata(hFig,'hTxt_Frame'));
ReSzFig(hFig);

hAx = getappdata(hFig,'hAx');
hCB = getappdata(hFig,'hCB');

UseCB = ishandle(hCB);
nFrames = getappdata(hFig,'nFrames');
Anim(nFrames) = struct('cdata',[],'colormap',[]);
for n=1:nFrames
    change_frame(hFig,n);
    drawnow;
    if ~UseCB
        Anim(n) = getframe(hAx);
    else
        hFig.Units = 'pixels';
        hAx.Units = 'pixels';
        hCB.Units = 'pixels';
        tiAx = hAx.TightInset;
        tiCB = hCB.TightInset;
        pos = [min(tiAx(1),tiCB(1)),...
            min(tiAx(2),tiCB(2)),...
            max(tiAx(1)+tiAx(3),tiCB(1)+tiCB(3)),...
            max(tiAx(2)+tiAx(4),tiCB(2)+tiCB(4))];
        pos(3) = pos(3)-pos(1);
        pos(4) = pos(4)-pos(2);
        pos = round(pos);
        Anim(n) = getframe(hFig);
    end
end

%% Save Movie
writerObj = VideoWriter(fullfile(mov_path,mov_file),'MPEG-4');
writerObj.FrameRate = p.Results.FrameRate;
writerObj.Quality = p.Results.Quality;
open(writerObj);
writeVideo(writerObj,Anim);
close(writerObj);

%% recreate controls
set(hFig,'units','characters');
fpos = get(hFig,'position');
hSld_Frame = uicontrol('Parent', hFig,'Style','slider',...
         'Units','characters',...
         'Position',[1, 0, fpos(3)-15-1, 2],...
         'Min',1,...
         'Max',nFrames,...
         'Value',1,...
         'SliderStep',[1/nFrames,10/nFrames],...
         'Tag','hSld_Frame',...
         'Callback',@SldFrame_callback);

hEdt_Frame = uicontrol('Parent',hFig,...
        'Style','edit',...
        'units','characters',...
        'Position',[fpos(3)-15+1,0.7,7,1.6],...
        'String','1',...
        'HorizontalAlignment','center',...
        'Callback',@EditFrame_callback);
hTxt_Frame = uicontrol('Parent',hFig,...
    'Style','text',...
    'units','characters',...
    'position',[fpos(3)-15+1+7+.5,0,5,2],...
    'HorizontalAlignment','left',...
    'string',sprintf('/%d',nFrames));
setappdata(hFig,'hSld_Frame',hSld_Frame);
setappdata(hFig,'hEdt_Frame',hEdt_Frame);
setappdata(hFig,'hTxt_Frame',hTxt_Frame);
ReSzFig(hFig);

PostSaveFcn = getappdata(hFig,'PostSaveFcn');
if ~isempty(PostSaveFcn)
    PostSaveFcn();
end

%% re-enable resize
set(hFig,'Resize',fRsz);
end