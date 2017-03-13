function [Y_CROP,X_CROP] = uicropstack(stack,varargin)

global y_crop;
global x_crop;

hFig = stackfig(stack,varargin{:});
title('Draw cropping window and close figure to accept.');

hrect = imrect2('Parent',get(hFig,'CurrentAxes'),...
    'Color','r',...
    'handlevisibility','callback',...
    'LimMode','manual');

set(hFig,'CloseRequestFcn',{@CloseFigReq,hrect});

waitfor(hFig);

Y_CROP = round(y_crop);
X_CROP = round(x_crop);
Y_CROP = max(1,min(size(stack,1),Y_CROP));
X_CROP = max(1,min(size(stack,2),X_CROP));



function CloseFigReq(hFig,~,hRect)
global y_crop;
global x_crop;
pos = get(hRect,'Position');
x_crop = [pos(1),pos(1)+pos(3)-1];
y_crop = [pos(2),pos(2)+pos(4)-1];
delete(hFig);


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



function hrect = imrect2(varargin)
%Create Dynamic rectangle on image (similar to imrect() included with
%image processing toolbox but built with the standard rectangle object)
% Input:
%   'Parent',hax - handle to parent axes (default = gca)
%   'Position',[x,y,w,h] - position of rectangle in axes units
%                          if not specified rbbox is called for user to
%                          select one graphically
%   'Color',ColorSpec - Color of rectangle edge and corner markers
%   'LineWidth',w - rectangle line width (default=0.5)
%   'LineStyle,LineSpec - rectangle line style
%   'MarkerSize',w - size of corner markers (default = 12*LineWidth)
%   'ResizeFcn',fcn or 'fcn' or {fcn,arg1...} - Function to call when
%               rectangle size is changed
%               The first argument of fcn will be hrect (the handle to the
%               rectange) use get(hrect,'position') to get position after
%               resize.
%   'LimMode','auto' or 'manual' - if 'manual' axis lims are fixed at
%               limits before imrect2 was called.
%   'HandleVisibility','on'(default)|'off'|'callback'
%               Visibility of object handle (see Rectangle properties)
%               This is useful if you want to have the rectangle persist
%               after something like plot(...) is called.
%               If the parent axes is set to 'NextPlot'='replacechildren'
%               then setting 'handlevisibility'='callback' will prevent
%               plot(...) from deleting the rectangle even if hold on is
%               not set.
%   'LockPosition',true|false(default) - specify if the position is locked
%               after the rectangle has been created. prevent user from
%               modifying position. 
%               LockPosition can be change after the rectangle is created
%               by changing the userdata.
%                   Example:
%                       ud = get(hrect,'userdata'); %get userdata
%                       ud.LockPosition = true; %change value
%                       set(hrect,'userdata',ud); %save userdata
% Output:
%   hrect - handle to the rectangle
%       Note: the hrect will contain userdata with the following elements
%           userdata.hPt_BL - handle to bottom left marker
%           userdata.hPt_BR - handle to bottom right marker
%           userdata.hPt_TL - handle to top left marker
%           userdata.hPt_TR - handle to top right marker
%           userdata.ResizeFcn - the resize function set above
%           userdata.LockPosition - (true/false) flag specifying if 
%                                   rectangle location is locked 
%        Also be aware that the ButtonDownFcn has been set for hrect and
%        the corner points (hPt...)
%==========================================================================
% Copyright 2015, Daniel Kovari. All rights reserved.

hrect = [];
%input parser
p = inputParser;
p.CaseSensitive = false;
addParameter(p,'Parent',[]);
addParameter(p,'Position',[]);
addParameter(p,'Color','k');
addParameter(p,'HandleVisibility','on', @(x) any(strcmpi(x,{'on','off','callback'})));
addParameter(p,'LineWidth',0.5);
addParameter(p,'LineStyle','-');
addParameter(p,'MarkerSize',[]);
addParameter(p,'ResizeFcn',[],@verify_fcn);
addParameter(p,'LimMode','auto',@(x) any(strcmpi(x,{'auto','manual'})));
addParameter(p,'LockPosition',false, @islogical);
parse(p,varargin{:});

hparent = p.Results.Parent;
position = p.Results.Position;
if ~isempty(hparent)
    if ~ishandle(hparent)||~strcmpi(get(hparent,'type'),'axes')
        error('hparent was not a valid axes handle');
    end
else
    hparent = gca;
end
if ~isempty(position)
    if ~isnumeric(position)||numel(position)~=4
        error('position must be [x,y,w,h]');
    end
end

hfig = get(hparent,'parent');

if isempty(position)
    axes(hparent);
    ptr = get(hfig,'pointer');
    set(hfig,'pointer','crosshair');
    switch(waitforbuttonpress)
        case 0 %user clicked the mouse
            r = rbbox;
            set(hfig,'pointer',ptr);
        case 1 %user pressed a key return without creating box
            set(hfig,'pointer',ptr);
            return;
    end
    %[x,y] = fig2axescoord([r(1);r(1)+r(3)],[r(2);r(2)+r(4)],hparent);
    %position = [x(1),y(1),x(2)-x(1),y(2)-y(1)];
    position = fig2axpos(r,hparent);
end


if strcmpi(p.Results.LimMode,'manual')
    xl = get(hparent,'xlim');
    yl = get(hparent,'ylim');
end
hrect = rectangle('parent',hparent,...
    'position',position,...
    'HandleVisibility',p.Results.HandleVisibility,...
    'LineWidth',p.Results.LineWidth,...
    'LineStyle',p.Results.LineStyle,...
    'EdgeColor',p.Results.Color);
if strcmpi(p.Results.LimMode,'manual')
    set(hparent,'xlim',xl);
    set(hparent,'ylim',yl);
end

MarkerSize = p.Results.MarkerSize;
if isempty(MarkerSize)
    MarkerSize = 12*p.Results.LineWidth;
end

%setup resize handle points on corners
userdata.hPt_BL = line(position(1),position(2),0,...
    'parent',hparent,...
    'HandleVisibility',p.Results.HandleVisibility,...
    'marker','s',...
    'MarkerSize',MarkerSize,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',p.Results.Color,...
    'ButtonDownFcn',{@BL_click,hrect});
userdata.hPt_BR = line(position(1)+position(3),position(2),0,...
    'parent',hparent,...
    'HandleVisibility',p.Results.HandleVisibility,...
    'marker','s',...
    'MarkerSize',MarkerSize,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',p.Results.Color,...
    'ButtonDownFcn',{@BR_click,hrect});
userdata.hPt_TL = line(position(1),position(2)+position(4),0,...
    'parent',hparent,...
    'HandleVisibility',p.Results.HandleVisibility,...
    'marker','s',...
    'MarkerSize',MarkerSize,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',p.Results.Color,...
    'ButtonDownFcn',{@TL_click,hrect});
userdata.hPt_TR = line(position(1)+position(3),position(2)+position(4),0,...
    'parent',hparent,...
    'HandleVisibility',p.Results.HandleVisibility,...
    'marker','s',...
    'MarkerSize',MarkerSize,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',p.Results.Color,...
    'ButtonDownFcn',{@TR_click,hrect});

userdata.ResizeFcn = p.Results.ResizeFcn;
userdata.LockPosition = p.Results.LockPosition;

set(hrect,'userdata',userdata);
set(hrect,'ButtonDownFcn',@drag_rect);
set(hrect,'DeleteFcn',@delete_rect);

function drag_rect(hrect,~)
ud = get(hrect,'userdata');
if ud.LockPosition %check if lock position
    return;
end
hax = get(hrect,'parent');
p = get(hrect,'position');
%dragrect doesn't work with normalized units
orig_units = get(get(hax,'parent'),'units');
if strcmpi(orig_units,'normalized');
    set(get(hax,'parent'),'units','pixels');
end
%drag rectangle
r = dragrect(ax2figpos(p,hax));
p = fig2axpos(r); %get position value in axes units
%return figure units to whatever they were
set(get(hax,'parent'),'units',orig_units);

%set positions
set(hrect,'position',p);
%corner points
set(ud.hPt_BL,'xdata',p(1),'ydata',p(2));
set(ud.hPt_BR,'xdata',p(1)+p(3),'ydata',p(2));
set(ud.hPt_TL,'xdata',p(1),'ydata',p(2)+p(4));
set(ud.hPt_TR,'xdata',p(1)+p(3),'ydata',p(2)+p(4));
ExecResizeFcn(hrect);

function [x,y] = fig2axescoord(x,y,hax)
%disp('f2ac')
if nargin<3
    hax = gca;
end
par_units = get(get(hax,'parent'),'units');
orig_units = get(hax,'units');
set(hax,'units',par_units);
pos = plotboxpos(hax);
xl = get(hax,'xlim');
yl = get(hax,'ylim');
x = (x-pos(1))/pos(3)*(xl(2)-xl(1))+xl(1);
y = (y-pos(2))/pos(4)*(yl(2)-yl(1))+yl(1);
set(hax,'units',orig_units);

function [x,y] = axes2figcoord(x,y,hax)
%disp('a2fc')
if nargin<3
    hax = gca;
end
par_units = get(get(hax,'parent'),'units');
orig_units = get(hax,'units');
set(hax,'units',par_units);
pos = plotboxpos(hax);
xl = get(hax,'xlim');
yl = get(hax,'ylim');
x = (x-xl(1))/(xl(2)-xl(1))*pos(3)+pos(1);
y = (y-yl(1))/(yl(2)-yl(1))*pos(4)+pos(2);
set(hax,'units',orig_units);

function p = ax2figpos(p,hax)
if nargin<2
    hax = gca;
end
[x,y] = axes2figcoord([p(1),p(1)+p(3)],[p(2),p(2)+p(4)],hax);
p = [x(1),y(1),x(2)-x(1),y(2)-y(1)];

function p = fig2axpos(p,hax)
if nargin<2
    hax = gca;
end
[x,y] = fig2axescoord([p(1);p(1)+p(3)],[p(2);p(2)+p(4)],hax);
p = [x(1),y(1),x(2)-x(1),y(2)-y(1)];

function BL_click(~,~,hrect)
ud = get(hrect,'userdata');
if ud.LockPosition %check if lock position
    return;
end
hax = get(hrect,'parent');
p = get(hrect,'position');
r = ax2figpos(p,hax);
r = rbbox(r,[r(1)+r(3),r(2)+r(4)]);
p = fig2axpos(r,hax);

%set positions
set(hrect,'position',p);
%corner points
set(ud.hPt_BL,'xdata',p(1),'ydata',p(2));
set(ud.hPt_BR,'xdata',p(1)+p(3),'ydata',p(2));
set(ud.hPt_TL,'xdata',p(1),'ydata',p(2)+p(4));
set(ud.hPt_TR,'xdata',p(1)+p(3),'ydata',p(2)+p(4));
ExecResizeFcn(hrect);

function BR_click(~,~,hrect)
ud = get(hrect,'userdata');
if ud.LockPosition %check if lock position
    return;
end
hax = get(hrect,'parent');
p = get(hrect,'position');
r = ax2figpos(p,hax);
r = rbbox(r,[r(1),r(2)+r(4)]);
p = fig2axpos(r,hax);
%set positions
set(hrect,'position',p);
%corner points
set(ud.hPt_BL,'xdata',p(1),'ydata',p(2));
set(ud.hPt_BR,'xdata',p(1)+p(3),'ydata',p(2));
set(ud.hPt_TL,'xdata',p(1),'ydata',p(2)+p(4));
set(ud.hPt_TR,'xdata',p(1)+p(3),'ydata',p(2)+p(4));
ExecResizeFcn(hrect);

function TL_click(~,~,hrect)
ud = get(hrect,'userdata');
if ud.LockPosition %check if lock position
    return;
end
hax = get(hrect,'parent');
p = get(hrect,'position');
r = ax2figpos(p,hax);
r = rbbox(r,[r(1)+r(3),r(2)]);
p = fig2axpos(r,hax);
%set positions
set(hrect,'position',p);
%corner points
set(ud.hPt_BL,'xdata',p(1),'ydata',p(2));
set(ud.hPt_BR,'xdata',p(1)+p(3),'ydata',p(2));
set(ud.hPt_TL,'xdata',p(1),'ydata',p(2)+p(4));
set(ud.hPt_TR,'xdata',p(1)+p(3),'ydata',p(2)+p(4));
ExecResizeFcn(hrect);

function TR_click(~,~,hrect)
ud = get(hrect,'userdata');
if ud.LockPosition %check if lock position
    return;
end
hax = get(hrect,'parent');
p = get(hrect,'position');
r = ax2figpos(p,hax);
r = rbbox(r,[r(1),r(2)]);
p = fig2axpos(r,hax);
%set positions
set(hrect,'position',p);
%corner points
set(ud.hPt_BL,'xdata',p(1),'ydata',p(2));
set(ud.hPt_BR,'xdata',p(1)+p(3),'ydata',p(2));
set(ud.hPt_TL,'xdata',p(1),'ydata',p(2)+p(4));
set(ud.hPt_TR,'xdata',p(1)+p(3),'ydata',p(2)+p(4));
ExecResizeFcn(hrect);

function ExecResizeFcn(hrect)
ud = get(hrect,'userdata');
if isempty(ud.ResizeFcn)
    return;
end

if ischar(ud.ResizeFcn)
    f = str2func(ud.ResizeFcn);
elseif iscell(ud.ResizeFcn)
    if ischar(ud.ResizeFcn{1})
        f = str2func(ud.ResizeFcn{1});
        f = @(x) f(x,ud.ResizeFcn{2:end});
    else
        f = @(x) ud.ResizeFcn{1}(x,ud.ResizeFcn{2:end});
    end
elseif isa(ud.ResizeFcn,'function_handle')
    f = ud.ResizeFcn;
else
    error('Something is wrong with ud.ResizeFcn');
end
f(hrect);

function stat = verify_fcn(f)
if isa(f,'function_handle')
    stat=true;
elseif ischar(f)
    stat = true;
elseif iscell(f)&&(isa(f{1},'function_handle')||ischar(f{1}))
    stat = true;
elseif isempty(f)
    stat = true;
else
    stat = false;
end

function delete_rect(hrect,~)
%disp('indelete')
ud = get(hrect,'userdata');
%delete corner points
try
    delete(ud.hPt_BL);
    delete(ud.hPt_BR);
    delete(ud.hPt_TL);
    delete(ud.hPt_TR);
catch
end
%delete the rectangle;
delete(hrect);

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
% Output variables:
%   pos:    four-element position vector, in same units as h
% Copyright 2010 Kelly Kearney

% Check input
if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    pos = axisPos; 
else
    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis
temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);