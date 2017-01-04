function TFMForceViewer(TFMdata,varargin)
%Traction Force Viewer
% Syntax:
%   If no arguments are specified the user is propmpted to open a saved
%   traction force data file (*.mat format).
%
%   Alternatively, you can specify a file or TFMdata struct and various
%   display formating parameters
%
% Inputs
%   TFMdata: a string specifying the file path to open OR TFMdata struct
%   generated by CalculateTFM_*(). If TFMdata=[] then user is prompted to
%   select a file.
%
% Parameters:
%   'MoviePath': Specify a path to save force map movie
%   'CloseAfterSave',true/false: if true (default), the viewer is closed
%                                after the movie is saved.
%   'CellImageCLim',[low,high]: image intensity color limits for cell
%                               images. Alternatively specify 'average' or
%                               'global' to automatically determine limits
%   'PlotStrain',false/true: plot vectors indicating gel strain
%   'StrainColor',rgb or color strings specifying default strain vect color
%   'SMAGLim',[low,high']/string: Specify color limits for |Stress| map
%   'SMAGColor',rgb: color overlay for |Stress| map

%% Parse Inputs
p = inputParser;
p.CaseSensitive = false;
addParameter(p,'MoviePath',[]);
addParameter(p,'CloseAfterSave',true);
addParameter(p,'CellImageCLim','average');
addParameter(p,'PlotStrain',false);
addParameter(p,'StrainColor',[255,215,0]/255,...
    @(x) ischar(x)&&...
    any(strcmpi(x,{'r','c','k','b','y','w','m','g','yellow','magenta','cyan','red','green','blue','white','black'}))...
    ||isnumeric(x)&&numel(x)==3&&all(x<=1)&&all(x>=0));
addParameter(p,'SMAGLim','global');
addParameter(p,'SMAGColor',[1,0,0],@(x) isnumeric(x)&&numel(x)==3&&all(x<=1)&&all(x>=0));
addParameter(p,'FigureSize',[],@(x) isempty(x)||isnumeric(x)&&numel(x)==2);
parse(p,varargin{:});

%% Load Data
persistent last_dir;

if nargin<1||isempty(TFMdata)
    [File,Dir] = uigetfile(fullfile(last_dir,'*.mat'),'Calculated TFM Data');
    if File==0
        return
    end
    if ~isempty(Dir)
        last_dir = Dir;
    end
    TFMdata = fullfile(Dir,File);
end

if ischar(TFMdata) %TFMdata is a string specifying a file
    TFMdata = load(TFMdata);
end

%% Validate TFMdata
fields = {'Time','Ostack','Vxx','Vyy','Vqx','Vqy','SF','dx','dy','E','v','PX_SCALE','SMAG'};
if any( ~isfield(TFMdata,fields))
    erfld = fields(~isfield(TFMdata,fields));
    error('TFMdata Missing Field: %s\n',erfld{:});
end

%number of image pixels corresponding to each point in SMAG
dW = TFMdata.dx/TFMdata.PX_SCALE;
dH = TFMdata.dy/TFMdata.PX_SCALE;

%Convert SMAG stack to anim struct
olAnim(size(TFMdata.SMAG,3)) = struct('AlphaData',[],'XData',[],'YData',[]);
[olAnim.XData] = deal([dW/2, dW*(size(TFMdata.SMAG,2)-1)-dW/2]);
[olAnim.YData] = deal([dH/2, dH*(size(TFMdata.SMAG,1)-1)-dH/2]);
for n=1:size(TFMdata.SMAG,3)
    olAnim(n).AlphaData = TFMdata.SMAG(:,:,n);
end

%% Create the overlay figure
[hFig,hAx,hCB] = overlay_animfig(TFMdata.Ostack,olAnim,...
    'OverlayColor',p.Results.SMAGColor,...
    'CLim',p.Results.CellImageCLim,...
    'ALim',p.Results.SMAGLim,...
    'colormap',gray(512),...
    'frameupdate_fn',@FrameChange);
ylabel(hCB,'|Stress| [Pa]');
set(hAx,'XColor','none','YColor','none');
set(hAx,'XTick',[],'YTick',[]);
%% Draw Quivers and pointers
hold(hAx,'on');
%hCnt = plot(hAx,TFMdata.cnt{1}(:,1), TFMdata.cnt{1}(:,2),'+k','markersize',10);
if p.Results.PlotStrain
    hQvr = quiver(reshape(TFMdata.Vxx+dW/2,[],1),...
                reshape(TFMdata.Vyy+dH/2,[],1),...
                reshape(TFMdata.Vqx(:,:,1),[],1),...
                reshape(TFMdata.Vqx(:,:,1),[],1),...
                0,'-','color',p.Results.StrainColor,...
                'LineWidth',1.5);
    legend(hAx,hQvr,'Strain');
end
%set axes limits to fit image, so that vectors dont cause scaling isues
xlim(hAx,[0,size(TFMdata.Ostack,2)]);
ylim(hAx,[0,size(TFMdata.Ostack,1)]);
%% Create TimeStamp & Scalebar
TIME_FONT_SIZE = 26;
yloc = 0.99; %top of text (frac of axis)
xloc = 0.01; %left of text (frac of axis)

%location and size of scalebar
%Scale bar size in px
SB_FONT_SIZE = 26;
SB_LENGTH = 20; %�m
SB_WIDTH = 6; %figure points
PX_LENGTH = SB_LENGTH/(TFMdata.PX_SCALE/10^-6);

SB_POS = [0.05,0.05]; %position of scalebar [x,y] relative to lower corner

YLIM = get(gca,'ylim');
XLIM = get(gca,'xlim');

switch get(gca,'xdir')
    case 'normal'
        xloc = XLIM(1)+xloc*(XLIM(2)-XLIM(1));
        SB_X = XLIM(1)+SB_POS(1)*(XLIM(2)-XLIM(1)) + [0,PX_LENGTH];
    case 'reverse'
        xloc = XLIM(2)-xloc*(XLIM(2)-XLIM(1));
        SB_X = XLIM(2)-SB_POS(1)*(XLIM(2)-XLIM(1)) + [0,-PX_LENGTH];
end
switch get(gca,'ydir')
    case 'normal'
        yloc = YLIM(1)+yloc*(YLIM(2)-YLIM(1));
        SB_Y = YLIM(1)+SB_POS(2)*(YLIM(2)-YLIM(1)) + [0,0];
    case 'reverse'
        yloc = YLIM(2)-yloc*(YLIM(2)-YLIM(1));
        SB_Y = YLIM(2)-SB_POS(2)*(YLIM(2)-YLIM(1)) + [0,0];
end

%Time Stamp
str = sprintf('Time: %04.01f min',0);
hTS = text(xloc,yloc,str,'Color','w','VerticalAlignment','top','HorizontalAlignment','Left','FontSize',TIME_FONT_SIZE);

%Plot ScaleBar
if strcmp(get(gca,'ydir'),'normal')
    text(mean(SB_X),mean(SB_Y)+2*SB_WIDTH,sprintf('%d �m',SB_LENGTH),'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
else
    text(mean(SB_X),mean(SB_Y)-2*SB_WIDTH,sprintf('%d �m',SB_LENGTH),'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
end
plot(SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);

%% Animation Function
function FrameChange(~,~,f)
    %set(hCnt,'XData',TFMdata.cnt{1}(:,1),...
    %         'YData',TFMdata.cnt{1}(:,2));
    if p.Results.PlotStrain
        set(hQvr,'UData',reshape(TFMdata.Vqx(:,:,f),[],1),...
                 'VData',reshape(TFMdata.Vqy(:,:,f),[],1));
    end
    
    set(hTS,'String',sprintf('Time: %04.01f min',(TFMdata.Time(f)-TFMdata.Time(1))/60));
end

%% Write movie if needed
if ~isempty(p.Results.MoviePath)
    WriteFn = getappdata(hFig,'ExecuteMovie_Fn');
    %make background white
    old_bg = hFig.Color;
    hFig.Color = 'w';
    if isempty(p.Results.FigureSize)
        set(hFig,'units','normalized','position',[0,0,1,1]);
    else
        set(hFig,'units','pixels','position',[0,0,p.Results.FigureSize(1),p.Results.FigureSize(2)]);
    end
    WriteFn(p.Results.MoviePath);
    if p.Results.CloseAfterSave
        delete(hFig);
        return;
    else
        hFig.Color = old_bg;
    end
end

end