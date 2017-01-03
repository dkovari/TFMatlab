function TFMForceViewer(TFMdata,varargin)

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'MoviePath',[]);
addParameter(p,'CloseAfterSave',true);
addParameter(p,'CellImageCLim','average');
addParameter(p,'PlotDisplacements',false);
addParameter(p,'SMAGLim','global');
parse(p,varargin{:});

%% Load Data
persistent last_dir;

if nargin<1
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
[olAnim.XData] = deal([dW/2, dW*size(TFMdata.SMAG,2)-dW/2]);
[olAnim.YData] = deal([dH/2, dH*size(TFMdata.SMAG,1)-dH/2]);
for n=1:size(TFMdata.SMAG,3)
    olAnim(n).AlphaData = TFMdata.SMAG(:,:,n);
end

%% Create the overlay figure
[hFig,hAx,hCB] = overlay_animfig(TFMdata.Ostack,olAnim,...
    'OverlayColor',[1,0,0],...
    'CLim',p.Results.CellImageCLim,...
    'ALim',p.Results.SMAGLim,...
    'colormap',gray(256),...
    'frameupdate_fn',@FrameChange);
ylabel(hCB,'|Stress| [Pa]');
%% Draw Quivers and pointers
hold(hAx,'on');
%hCnt = plot(hAx,TFMdata.cnt{1}(:,1), TFMdata.cnt{1}(:,2),'+k','markersize',10);
if p.Results.PlotDisplacements
    hQvr = quiver(reshape(TFMdata.Vxx+dW/2,[],1),...
                reshape(TFMdata.Vyy+dH/2,[],1),...
                reshape(TFMdata.Vqx(:,:,1),[],1),...
                reshape(TFMdata.Vqx(:,:,1),[],1),...
                0,'-y');
end
%% Create TimeStamp & Scalebar
TIME_FONT_SIZE = 26;
yloc = 0.99; %top of text (frac of axis)
xloc = 0.01; %left of text (frac of axis)

%location and size of scalebar
%Scale bar size in px
SB_FONT_SIZE = 26;
SB_LENGTH = 20; %µm
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
    text(mean(SB_X),mean(SB_Y)+2*SB_WIDTH,sprintf('%d µm',SB_LENGTH),'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
else
    text(mean(SB_X),mean(SB_Y)-2*SB_WIDTH,sprintf('%d µm',SB_LENGTH),'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
end
plot(SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);

%% Animation Function
function FrameChange(~,~,f)
    %set(hCnt,'XData',TFMdata.cnt{1}(:,1),...
    %         'YData',TFMdata.cnt{1}(:,2));
    if p.Results.PlotDisplacements
        set(hQvr,'UData',reshape(TFMdata.Vqx(:,:,f),[],1),...
                 'VData',reshape(TFMdata.Vqy(:,:,f),[],1));
    end
    
    set(hTS,'String',sprintf('Time: %04.01f min',(TFMdata.Time(f)-TFMdata.Time(1))/60));
end

%% Write movie if needed
if ~isempty(p.Results.MoviePath)
    WriteFn = getappdata(hFig,'ExecuteMovie_Fn');
    set(hFig,'units','normalized','position',[0,0,1,1]);
    WriteFn(p.Results.MoviePath);
    if p.Results.CloseAfterSave
        delete(hFig);
        return;
    end
end

end