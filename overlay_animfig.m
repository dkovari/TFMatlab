function [hFig,hAx,hCB] = overlay_animfig(Base,Over,varargin)
% Overlay an AlphaMap "movie" on top of a ColorMapped "movie"
% This view utilizes composit_anigfig(...)
%
% Input:
%   Base: Color or grayscale image sequence ontop of which AlphaMap will be
%         displayed. If image is grayscale, it will converted to an RGB
%         image using the specified colormap
%
%   Over: Grayscale image sequence to overlay.
%
% Both Base and Over can be specified as numeric arrays, cell arrays or
% animation structures (see composit_animfig)
%
% The color of Over can be specified by 'OverlayColor', or by specifying an
% image in the anim struct Over(n).CData = ...
%
% Parameters:
%   'colormap': colormap to use for Base images
%   'CLim': colorlims to use when applying colormap (see composit_animfig)
%   'ALim': Alim of axes
%       'auto': use the default values (scales to match image extents)
%       'global': use max/min for all frames in Over
%       'average' use average of max/min for frames in Over
%   'Figure',hFig: specify figure to use;
%   'frameupdate_fn':@(hFig,hAx,currentFrame) specify a callback to execute when frame
%                    is changed. See compsit_animfig for how to retrieve
%                    current frame and image data.

if ~isstruct(Over) && ~iscell(Over) && ~(isnumeric(Over)&&ndims(Over)<4)
    error('Invalid type for over');
end

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'colormap','gray',@(x) ischar(x)||isnumeric(x)&&(size(x,2)==3));
addParameter(p,'OverlayColor',[1,0,0], @(x) isnumeric(x)&&numel(x)==3);
addParameter(p,'ALim','global',@(x) ischar(x)&&any(strcmpi(x,{'auto','global','average'}))||isnumeric(x)&&numel(x)==2);
addParameter(p,'CLim','global',@(x) ischar(x)&&any(strcmpi(x,{'direct','scaled','global','average'}))||isnumeric(x)&&numel(x)==2);
addParameter(p,'Figure',[],@(x) isempty(x)||ishghandle(x));
addParameter(p,'frameupdate_fn',[],@(x) isempty(x)||isa(x,'function_handle'));

parse(p,varargin{:});

if iscell(Over)
    cOver = Over;
    Over = struct('AlphaData',{});
    for n=1:numel(cOver)
        if ~isnumeric(cOver{n})||ndims(cOver{n})>2
            error('Specified cell array for Over contained non-image (matrix) data');
        end
        Over(n).AlphaData = cOver;
    end
end

if isnumeric(Over)
    mOver = Over;
    Over = struct('AlphaData',{});
    for n=1:size(mOver,3)
        Over(n).AlphaData = mOver(:,:,n);
    end
end

%% Data is struct, add color information
if ~isfield(Over,'AlphaData')
    error('If specifying Over as a struct, it must contain AlphaData')
end
if ~isfield(Over,'CData')
    for n=1:numel(Over)
        [H,W] = size(Over(n).AlphaData);
        Over(n).CData = repmat(reshape(p.Results.OverlayColor,1,1,[]),H,W);
    end
end
%% Determine ALim
ALim = p.Results.ALim;
Alow = NaN(numel(Over),1);
Ahigh = NaN(numel(Over),1);
for n=1:numel(Over)
    Alow(n) = nanmin(Over(n).AlphaData(:));
    Ahigh(n) = nanmax(Over(n).AlphaData(:));
end

if ischar(ALim)
    switch lower(ALim)
        case'auto'

        case 'global'
            ALim = [min(Alow),max(Ahigh)];
        case 'average'
            ALim = [nanmean(Alow),nanmean(Ahigh)];
    end
end

AExt = [min(Alow),max(Ahigh)];

% All overdata should be scaled to the alphamap
[Over.AlphaDataMapping] = deal('scaled');

%% Create AnimFig
[hFig,hAx,~,hCB] = composit_animfig(Base,Over,...
                            'colormap',p.Results.colormap,...
                            'CLim',p.Results.CLim,...
                            'Figure',p.Results.Figure,...
                            'frameupdate_fn',@FrameUpdateFcn,...
                            'ShowColorbar',true,...
                            'ColorbarWidth',20,...
                            'PreSaveFcn',@PreSaveAnim,...
                            'PostSaveFcn',@PostSaveAnim);
%% Setup ALim
if isnumeric(ALim)
    set(hAx,'ALim',ALim);
    set(hCB,'ALim',ALim);
else
    set(hAx,'ALimMode','auto');
    set(hCB,'ALimMode','auto');
end
                        
%% create colorscale image & histogram in the colorbar

[Counts,edges] = histcounts(Over(1).AlphaData(:),'BinLimits',AExt);
[Y,X] = edges2stairs(edges,Counts);
lX = log10(X);
lX(X==0) = 0;
hHistLine = plot(hCB,lX,Y,'-k','linewidth',1.5,'hittest','off');
hold(hCB,'on');

XLIM = get(hCB,'XLim');
ylim(hCB,AExt);
cb_cdata = repmat(reshape(p.Results.OverlayColor,1,1,[]),200,1);
cb_adata = linspace(AExt(1),AExt(2),200)';
hImCB = image(hCB,'CData',cb_cdata,...
    'AlphaData',cb_adata,...
    'AlphaDataMapping','scaled',...
    'XData',[0,XLIM(2)],...
    'YData',AExt,...
    'HitTest','off');

axis(hCB,'tight');
set(hCB,'XLim',[0,XLIM(2)]);
uistack(hImCB,'down');

AL = get(hAx,'ALim');
hL_low = plot(hCB,[0,XLIM(2)],[AL(1),AL(1)],':k','linewidth',3);
hL_up = plot(hCB,[0,XLIM(2)],[AL(2),AL(2)],':k','linewidth',3);

set(hCB,'ALim',AL);

set(hL_low,'ButtonDownFcn',@(h,~) LowBtnDwn(h,hCB,hAx,hImCB));
set(hL_up,'ButtonDownFcn',@(h,~) UpBtnDwn(h,hCB,hAx,hImCB));

set(hCB,'XDir','Reverse',...
    'XTick',[],...
    'Box','off',...
    'YAxislocation','right',...
    'TickDir','out');
xlabel(hCB,'Log_{10}(Count)');

    function FrameUpdateFcn(hF,~,cF)
        %cF = getappdata(hF,'curFrame');
        FR = getappdata(hF,'FrameData');
        od = FR{2}(cF).AlphaData;
        [C,E] = histcounts(od(:),'BinLimits',AExt);
        [Y,X] = edges2stairs(E,C);
        lX = log10(X);
        lX(X==0) = 0;
        set(hHistLine,'xdata',lX,'ydata',Y);
        
        if ~isempty(p.Results.frameupdate_fn)
            p.Results.frameupdate_fn(hF,getappdata(hF,'hAx'),cF)
        end
    end

    function PreSaveAnim()
        ylim(hCB,get(hAx,'ALim'));
        set(hHistLine,'visible','off');
        set(hL_low,'visible','off');
        set(hL_up,'visible','off');
        xlabel(hCB,'Log_{10}(Count)');
        delete(hCB.XLabel);
    end
    function PostSaveAnim()
        ylim(hCB,AExt);
        set(hHistLine,'visible','on');
        set(hL_low,'visible','on');
        set(hL_up,'visible','on');
        xlabel(hCB,'Log_{10}(Count)');
    end
end

function LowBtnDwn(h,hCB,hAx,hImCB)
set(gcf,'WindowButtonUpFcn',@RelFn);
set(gcf,'WindowButtonMotionFcn',@MotFn)
    function MotFn(~,~)
        cp = hCB.CurrentPoint;
        y = cp(1,2);
        set(h,'ydata',[y,y]);
    end
    function RelFn(hF,~)
        set(hF,'WindowButtonUpFcn',[]);
        set(hF,'WindowButtonMotionFcn',[]);
        AL = get(hAx,'ALim');
        AL(1) = h.YData(1);
        set(hAx,'ALim',AL);
        set(hCB,'ALim',AL);
        Ext = get(hCB,'YLim');
        cb_adata = linspace(Ext(1),Ext(2),200)';
        set(hImCB,'YData',Ext,'AlphaData',cb_adata);
        
    end

end

function UpBtnDwn(h,hCB,hAx,hImCB)
set(gcf,'WindowButtonUpFcn',@RelFn);
set(gcf,'WindowButtonMotionFcn',@MotFn)
    function MotFn(~,~)
        cp = hCB.CurrentPoint;
        y = cp(1,2);
        set(h,'ydata',[y,y]);
    end
    function RelFn(hF,~)
        set(hF,'WindowButtonUpFcn',[]);
        set(hF,'WindowButtonMotionFcn',[]);
        AL = get(hAx,'ALim');
        AL(2) = h.YData(1);
        set(hAx,'ALim',AL);
        set(hCB,'ALim',AL);
        Ext = get(hCB,'YLim');
        cb_adata = linspace(Ext(1),Ext(2),200)';
        set(hImCB,'YData',Ext,'AlphaData',cb_adata);
    end

end