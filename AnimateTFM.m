function AnimateTFM(TFMdata,MovDir,MovName)
persistent last_dir;
if nargin<1
    %% Prompt User for data
    %select file
    [File,Dir] = uigetfile(fullfile(last_dir,'*.mat'),'Select .mat file with TFM data');
    if File==0
        return
    end
    if ~isempty(Dir)
        last_dir = Dir;
    end
    TFMdata = load(fullfile(Dir,File));
    [~,name,~] = fileparts(File);
else
    if ischar(TFMdata)
        if ~exist(TFMdata,'file')
            error('Could not find the specified file');
        end
        [Dir,name,~] = fileparts(TFMdata);
        TFMdata = load(TFMdata);
    else
        if ~isstruct(TFMdata)
            error('TFMdata must be a struct containing tfm variables, or a string specifying the location of a file containing those variables');
        end
        Dir = '';
        name = '*';
    end
end

%% Validate file

if ~isfield(TFMdata,'PX_SCALE')
    error('Data does not include PX_SCALE');
end
if ~isfield(TFMdata,'Ostack')
    error('Data does not include Ostack');
end
if ~isfield(TFMdata,'SMAG')
    error('Data does not include SMAG');
end
if ~isfield(TFMdata,'bpass_lnoise')
    error('Data does not include bpass_lnoise');
end
if ~isfield(TFMdata,'bpass_sz')
    error('Data does not include bpass_sz');
end
if ~isfield(TFMdata,'imstack')
    error('Data does not include imstack');
end
if ~isfield(TFMdata,'tracks')
    error('Data does not include tracks');
end
if ~isfield(TFMdata,'disptracks')
    error('Data does not include disptracks');
end
if ~isfield(TFMdata,'Vxx')
    error('Data does not include Vxx');
end
if ~isfield(TFMdata,'Vyy')
    error('Data does not include Vqx');
end
if ~isfield(TFMdata,'Vyy')
    error('Data does not include Vqx');
end
if ~isfield(TFMdata,'Vqy')
    error('Data does not include Vqy');
end
if ~isfield(TFMdata,'SF')
    error('Data does not include SF');
end
if ~isfield(TFMdata,'Time')
    error('Data does not include Time');
end


%% Other Parameters
PLOT_PART_DISP = false;
PLOT_DISP_FIELD = true;
PLOT_FORCE_VEC = false;
SHOW_BEADS = false;

if ~isfield(TFMdata,'cnt')&&PLOT_PART_DISP
    error('Data does not include cnt');
end

%% Initialize Figure
[H,W,nF] = size(TFMdata.Ostack);
PX_SCALE = TFMdata.PX_SCALE;
FSCALE = 10/100; %px/Pa
DSCALE = 2;
clim=stackclim(TFMdata.Ostack,'average');
cmap = gray(255);
alphaSMAG = 1;
climSMAG = [0,1500];%stackclim(TFMdata.SMAG,'global');
colorSMAG = [1,0,0];

%Plot Time Stamp
SHOW_TIMESTAMP = true;
%Time Stamp Placement
TIME_FONT_SIZE = 28;
yloc = 0.99; %top of text (frac of axis)
xloc = 0.01; %left of text (frac of axis)

%location and size of scalebar
%Scale bar size in px
SB_FONT_SIZE = 28;
SB_LENGTH = 20; %µm
SB_WIDTH = 6; %figure points
PX_LENGTH = SB_LENGTH/(PX_SCALE/10^-6);

RELATIVE_SB  = true; %place scalebar relative to axes corner
SB_POS = [0.05,0.05]; %position of scalebar [x,y]
SB_X = [10,10 + PX_LENGTH]; %non-relative position of SB
SB_Y = [30,30];%non-relative position of SB



%% setup fig and axes
hfig = figure('units','pixels','Position',[0,0,2*W+20,2*H+20]);
hax = axes('Parent',hfig);

%setup image
imbase = ind2rgb( gray2ind( mat2gray(TFMdata.Ostack(:,:,1),clim),size(cmap,1)),cmap);
hImage = image('Parent',hax,'CData',imbase,'handlevisibility','off');
axis(hax,'xy','image');

set(hax,'box','off',...
    'xtick',[],...
    'ytick',[]);

if RELATIVE_SB
    YLIM = get(hax,'ylim');
    XLIM = get(hax,'xlim');
    SB_X = XLIM(1)+SB_POS(1)*(XLIM(2)-XLIM(1)) + [0,PX_LENGTH];
    SB_Y = YLIM(1)+SB_POS(2)*(YLIM(2)-YLIM(1)) + [0,0];
end
 
%generate colormap
meanRGB(1) = mean2(imbase(:,:,1));
meanRGB(2) = mean2(imbase(:,:,2));
meanRGB(3) = mean2(imbase(:,:,3));
mnRGB = cat(3,meanRGB(1)*ones(100,1),meanRGB(2)*ones(100,1),meanRGB(3)*ones(100,1));
lsSMAG = linspace(climSMAG(1),climSMAG(2),100)';
%create colormap from overlay values
%CBcmap = imoverlaycmap(mnRGB,lsSMAG,cmapSMAG,alphaSMAG,climSMAG);
CBcmap = imOverlayOnRGB(mnRGB,lsSMAG,colorSMAG,climSMAG,[0,alphaSMAG]);
CBcmap = permute(CBcmap,[1,3,2]);

%color lims
set(hax,'CLim',climSMAG);
%set(hax,'CLimMode','manual');
%set colormap
%CBcmap=jet(100);
colormap(hax,CBcmap);

%% show colorbar
hcb = colorbar(hax,'location','East');
%colorbar ticks
nTICKS = 9;
CBticks = linspace(climSMAG(1),climSMAG(2),nTICKS);
set(hcb,'Ticks',CBticks);
% set(hcb,'YLimMode','manual');
% set(hcb,'YLim',[1,100]);
% set(hcb,'YTickMode','manual');
% set(hcb,'YTick',linspace(1,100,nTICKS));
% set(hcb,'YTickLabel',num2str(CBticks','%0.1f'));

%set colorbar font
set(hcb,'FontWeight','bold');
set(hcb,'FontSize',16);
set(hcb,'YColor','w');
set(hcb,'XColor','w');
%colorbar title
hcbtitle = get(hcb,'ylabel');
set(hcbtitle,'string','Stress Magnitude [Pa]');
set(hcbtitle,'FontSize',16);
set(hcbtitle,'FontWeight','Bold');


% Time stamp
YLIM = get(gca,'ylim');
XLIM = get(gca,'xlim');

switch get(gca,'xdir')
    case 'normal'
        xloc = XLIM(1)+xloc*(XLIM(2)-XLIM(1));
    case 'reverse'
        xloc = XLIM(2)-xloc*(XLIM(2)-XLIM(1));
end
switch get(gca,'ydir')
    case 'normal'
        yloc = YLIM(1)+yloc*(YLIM(2)-YLIM(1));
    case 'reverse'
        yloc = yLIM(2)-yloc*(YLIM(2)-YLIM(1));
end

Anim(nF) = struct('cdata',[],'colormap',[]);
%pause;
for f=1:nF
    cla(hax);
    
    %make rgb data
    imbase = ind2rgb( gray2ind( mat2gray(TFMdata.Ostack(:,:,f),clim),size(cmap,1)),cmap);
    
    %bead image
    if SHOW_BEADS
        B = bpass(TFMdata.imstack(:,:,f),TFMdata.bpass_lnoise,TFMdata.bpass_sz); 
        imbase = imOverlayOnRGB(imbase,B,[0,1,1],[0,max(B(:))*0.75],[0,.75]);
    end
    
    imSMAG = imresize(TFMdata.SMAG(:,:,f),[size(imbase,1),size(imbase,2)],'nearest');
    %im = imoverlaycmap(imbase,imSMAG,cmapSMAG,SMAGalpha,climSMAG);
    im = imOverlayOnRGB(imbase,imSMAG,colorSMAG,climSMAG,[0,alphaSMAG]);
    
    %show rgb image
    set(hImage,'CData',im);
    
    axis(hax,'xy','image');
    hold(hax,'on');
    
    %plot particle displacements
    if PLOT_PART_DISP
        plot(hax,reshape(TFMdata.tracks(f,1,:),[],1),reshape(TFMdata.tracks(f,2,:),[],1),'+k','markersize',10);
        plot(hax,TFMdata.cnt{f}(:,1), TFMdata.cnt{f}(:,2),'xb','markersize',8);

        quiver(hax,reshape(TFMdata.tracks(TFMdata.rF,1,:),[],1),...
                reshape(TFMdata.tracks(TFMdata.rF,2,:),[],1),...
                reshape(TFMdata.disptracks(f,1,:),[],1),...
                reshape(TFMdata.disptracks(f,2,:),[],1),...
                0,'-m');
    end
    %plot(hax,VxxPTV(isnan(Vqxn(:,:,f))),VyyPTV(isnan(Vqxn(:,:,f))),'ow');
    
    %plot displacement field
    if PLOT_DISP_FIELD
        quiver(hax,TFMdata.Vxx,TFMdata.Vyy,...
        DSCALE*TFMdata.Vqx(:,:,f),...
        DSCALE*TFMdata.Vqy(:,:,f),...
        0,'-y');
    end
    
    %plot calculated force
    if PLOT_FORCE_VEC
        quiver(hax,TFMdata.Vxx,TFMdata.Vyy,FSCALE*TFMdata.SF(:,:,1,f),FSCALE*TFMdata.SF(:,:,2,f),0,'-y');
    end
    
    %Time Stamp
    if SHOW_TIMESTAMP
        str = sprintf('Time: %04.01f min',(TFMdata.Time(f)-TFMdata.Time(1))/60);
        text(xloc,yloc,str,'parent',hax,'Color','w','VerticalAlignment','top','HorizontalAlignment','Left','FontSize',TIME_FONT_SIZE);
    end
    
    %Plot ScaleBar
    text(mean(SB_X),mean(SB_Y)+2*SB_WIDTH,sprintf('%d µm',SB_LENGTH),'parent',hax,'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
    plot(hax,SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);
    
    %save frame to animation
    xlim(hax,XLIM);
    ylim(hax,YLIM);
    
    Anim(f) = getframe(hax);
end
try
close(hfig);
catch
end

%putvar(Anim);
%implay(Anim);

%% Save mp4
if nargin>1
    if nargin>2
        [~,name,~] = fileparts(MovName);
    else
        if strcmp(name,'*')
            name = 'TFM';
        end
    end
    mov_file =[name,'.mp4'];
    mov_path = MovDir;
else
    [mov_file, mov_path] = uiputfile(fullfile(Dir,[name,'.mp4']),'Save Animation?');
end
if mov_file~=0
    %crop anim data to same size (this is a hack)
    Anim = SameSizeAnim(Anim);
    
    writerObj = VideoWriter(fullfile(mov_path,mov_file),'MPEG-4');
    writerObj.FrameRate = 5;
    open(writerObj);
    writeVideo(writerObj,Anim);
    close(writerObj);
end


