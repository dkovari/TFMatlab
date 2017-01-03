function hFig = stackfig(stack,varargin)
% Show stack on a figure with a scroll bar which can be used to cycle
% through images
% ========================================================================
% Inputs:
%   stack [H,W,nFrames] or [H,W,RGB,nFrames]
% Optional Arguments:
%   'CLim'
%       'auto' autoscales image mapping to [min,max] of each frame (default)
%       'global' uses [min,max] for all images
%       'average' calculates an average [min,max] for the set
%       [min,max] set clim to a specific value
%       clim[nF,2]  sets clim for each frame
%   'Colormap' specify the colormap using same syntax as colormap() in
%           matlab help
%   'FigureHandle', hFig
%           Specify the figure handle
%   'frameupdate_fn',@myfn() provide a function handle that will be called
%                           when ever the user changes the frame
%                           function should have the syntax
%                           myfn(hObject)
%                           where hObject is the handle to the figure and
%                           contains the data for the gui:
%                           handles = guidata(hObject); (see matlab help)
%           handles is a structure with the field
%             handles.curFrame: the current frame
%             handles.nFrames: the number of frames
%             handles.stack: the stack data
%             handles.RGBSTACK: bool specifying if stack is RGB or not
%             handles.AUTOSCALE: bool if using autoscale for image levels
%             handles.hImg: handle to the image being displayed
%             handles.hAx_Img: handle to the axis containg image
%             handles.hSld_Frame: handle to the slider
%             handles.hFig: handle to the figure
%             handles.userdata: userdata passed by the 'userdata' option
%             handles.frameupdate_fn: user specified frame update callback
%             handles.CLIM: clim to use for images
%             handles.CLIM_PERFRAME: bool specifying if CLIM is for single frames
%             handles.hTxt_Frames: handle to text specifying frame number
%                       
%                           
% Output:
%   hFig = handle to the figure
% *************************************************************************
% Copyright 2014 by Daniel T. Kovari
% Georgia Institute of Technology, School of Physics
% All rights reserved.

%% Check that stack is correct size
if ismatrix(stack)
    error('stack must have 3 dimensions');
elseif ndims(stack)==3
    RGBSTACK = false;
elseif ndims(stack)>4
    error('stack has too many dimensions');
else
    if size(stack,3)~=3
        error('RGB stack must have 3 color channels in dim 3');
    end
    RGBSTACK = true;
end

if RGBSTACK
    nFrames = size(stack,4);
else
    nFrames = size(stack,3);
end

%% parse variable inputs
p = inputParser;
p.CaseSensitive = false;

addParameter(p,'userdata',[]);
addParameter(p,'CLim','auto');
addParameter(p,'frameupdate_fn',[],@(x) isempty(x)||isa(x,'function_handle'));
addParameter(p,'Colormap',[],@(x) isnumeric(x)||ischar(x));
addParameter(p,'Figure',[],@(x) isempty(x)||ishghandle(x));

parse(p,varargin{:});

% validate CLim
CLIM = p.Results.CLim;
CLIM_PERFRAME = false;
AUTOSCALE = true;
if ~RGBSTACK  %don't bother with clim if stack is in color
    if isnumeric(CLIM)
        if length(CLIM)==nFrames
            CLIM_PERFRAME = true;
            if size(CLIM,2)~=2
                error('CLim must be either [nFrames,2] in size or have 2 elements');
            end
        else
            CLIM_PERFRAME = false;
            AUTOSCALE = false;
            if numel(CLIM)~=2
                warn('CLim must be either [nFrames,2] in size or have 2 elements, using auto instead');
                AUTOSCALE = true;
            end
        end
    elseif ischar(CLIM)
        CLIM = lower(CLIM);
        switch(CLIM)
            case 'auto'
                AUTOSCALE = true;
            case 'global'  %calculate global limits
                AUTOSCALE = false;
                CLIM = [nanmin(stack(:)), nanmax(stack(:))];
            case 'average'  %calc average limits
                AUTOSCALE = false;
                Maxlist = nan(nFrames,1);
                Minlist = Maxlist;
                for f = 1:nFrames
                    Maxlist(f) = nanmax(reshape(stack(:,:,f),[],1));
                    Minlist(f) = nanmin(reshape(stack(:,:,f),[],1));
                end
                CLIM = [nanmean(Minlist),nanmean(Maxlist)];
            otherwise
                warn('invalid clim method, using auto instead');
                AUTOSCALE = TRUE;
        end
    end
else
    AUTOSCALE = true;
end

%validate colormap
cmap = p.Results.Colormap;
USECMAP = false;
if isempty(cmap)
    USECMAP = false;
else
    USECMAP = true;
    if isnumeric(cmap)
        if size(cmap,2)~=3
            warn('colormap not specified correctly, see colormap() in matlab help.');
            USECMAP = false;
        end
    end
end          


%% initialize figure
if isempty(p.Results.Figure)
    hFig = figure();
    clf(hFig,'reset');
else
    clf(p.Results.Figure,'reset');
    hFig = figure(p.Results.Figure);
end

hSld_Frame = uicontrol('Parent', hFig,'Style','slider',...
         'Units','normalized',...
         'Position',[0, 0, 1, 0.05],...
         'Min',1,...
         'Max',nFrames,...
         'Value',1,...
         'SliderStep',[1/nFrames,10/nFrames],...
         'Tag','hSld_Frame',...
         'Callback',@SldFrame_callback);
%draw frame counter
hTxt_Frame = uicontrol('Style','text',...
                        'Units','normalized',...
                        'Position',[.9,.95,.1,.05],...
                        'String',['1','/',num2str(nFrames)],...
                        'BackgroundColor',get(hFig,'Color'));

hAx_Img = axes('Parent',hFig,...
           'Tag','hAx_Img',...
           'OuterPosition',[0,0.05,1,.95]);
if ~RGBSTACK
    hImg = imagesc('Parent',hAx_Img,'CData',stack(:,:,1));
    if ~AUTOSCALE
        if CLIM_PERFRAME
            set(hAx_Img,'CLim',CLIM(1,:));
        else
            set(hAx_Img,'CLim',CLIM);
        end
    end
else
    hImg = image('Parent',hAx_Img,'CData',stack(:,:,:,1));
end
axis(hAx_Img,'image');
if USECMAP
    colormap(hAx_Img,cmap);
end
set(hFig,'Toolbar','figure');



%% setup the guidata
handles = guidata(hFig);
handles.curFrame = 1;
handles.nFrames = nFrames;
handles.stack = stack;
handles.RGBSTACK = RGBSTACK;
handles.AUTOSCALE = AUTOSCALE;
handles.hImg =hImg;
handles.hAx_Img = hAx_Img;
handles.hSld_Frame = hSld_Frame;
handles.hFig = hFig;
handles.userdata = p.Results.userdata;
handles.frameupdate_fn = p.Results.frameupdate_fn;
handles.CLIM = CLIM;
handles.CLIM_PERFRAME = CLIM_PERFRAME;
handles.hTxt_Frames = hTxt_Frame;
guidata(hFig,handles);

end
%% Callback functions
function SldFrame_callback(hObject,~)
    %get new frame value
    handles = guidata(hObject);
    handles.curFrame = round(get(hObject,'Value'));
    %update cdata in image
    if ~handles.RGBSTACK %adjust limits for grayscale image
        set(handles.hImg,'CData',handles.stack(:,:,handles.curFrame));
        if handles.AUTOSCALE
            im = handles.stack(:,:,handles.curFrame);
            lim = [nanmin(im(:)),nanmax(im(:))];
            if lim(1)==lim(2)
                lim = [lim(1),lim(1)+1];
            end
            set(handles.hAx_Img,'CLim',lim);
        else
            if handles.CLIM_PERFRAME
                set(handles.hAx_Img,'CLim',handles.CLIM(handles.curFrame,:));
            else
                set(handles.hAx_Img,'CLim',handles.CLIM);
            end
        end
    else
        set(handles.hImg,'CData',handles.stack(:,:,:,handles.curFrame));
    end
    set(handles.hTxt_Frames,'String',[num2str(handles.curFrame),'/',num2str(handles.nFrames)]);

guidata(hObject,handles);

    if(~isempty(handles.frameupdate_fn))
        handles.frameupdate_fn(hObject);
    end
end