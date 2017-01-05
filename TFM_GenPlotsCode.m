%% Generate Plots of Results
[mov_file, mov_path] = uiputfile(fullfile(Dir,[name,sprintf('_XY%02d.mp4',SeriesNum)]),'Save Animation?');
if mov_file==0
    return;
end

PX_SCALE = 0.107498 * 1e-6; %um/px
FSCALE = 10/100; %px/Pa
DSCALE = 1;
clim=stackclim(Ostack,'average');
cmap = gray(255);
alphaSMAG = 1;
climSMAG = stackclim(SMAG,'global');
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

[H,W,nF] = size(Ostack);

%% setup fig and axes
hfig = figure('Position',[0,0,2*W+20,2*H+20]);
hax = axes('Parent',hfig);

%setup image
imbase = ind2rgb( gray2ind( mat2gray(Ostack(:,:,1),clim),size(cmap,1)),cmap);
image('Parent',hax,'CData',imbase);
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
set(hax,'CLimMode','manual');
%set colormap
%CBcmap=jet(100);
colormap(hax,CBcmap);

%% show colorbar
hcb = colorbar('peer',hax,'location','East');
%colorbar ticks
nTICKS = 9;
CBticks = linspace(climSMAG(1),climSMAG(2),nTICKS);
set(hcb,'YLimMode','manual');
set(hcb,'YLim',[1,100]);
set(hcb,'YTickMode','manual');
set(hcb,'YTick',linspace(1,100,nTICKS));
set(hcb,'YTickLabel',num2str(CBticks','%0.1f'));

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

if exist('Anim','var')
    clear Anim
end
Anim(nF) = struct('cdata',[],'colormap',[]);
%pause;
for f=1:nF
    cla(hax);
    
    %make rgb data
    imbase = ind2rgb( gray2ind( mat2gray(Ostack(:,:,f),clim),size(cmap,1)),cmap);
    
    %bead image
    B = bpass(imstack(:,:,f),bpass_lnoise,bpass_sz); 
    imbase = imOverlayOnRGB(imbase,B,[0,1,1],[0,max(B(:))*0.75],[0,.75]);
    
    imSMAG = imresize(SMAG(:,:,f),[size(imbase,1),size(imbase,2)],'nearest');
    %im = imoverlaycmap(imbase,imSMAG,cmapSMAG,SMAGalpha,climSMAG);
    im = imOverlayOnRGB(imbase,imSMAG,colorSMAG,climSMAG,[0,alphaSMAG]);
    
    
    
    %im=imbase;
    %im=ind2rgb( gray2ind( mat2gray(imSMAG,climSMAG),256),gray(256));
    %show rgb image
    image('Parent',hax,'CData',im);

    %imagesc('Parent',hax,'CData',Ostack(:,:,f));
    %colormap(hax,cmap);
    %set(hax,'CLim',clim);
    axis(hax,'xy','image');
    hold(hax,'on');
%    %% plot particle displacements
%     plot(reshape(tracks(1,1,:),[],1),reshape(tracks(1,2,:),[],1),'+m');
    %plot(reshape(tracks(f+1,1,:),[],1),reshape(tracks(f+1,2,:),[],1),'+k','markersize',10);
    
    plot(cnt{f+1}(:,1), cnt{f+1}(:,2),'+k','markersize',10);
    
    quiver(reshape(tracks(1,1,:),[],1),...
            reshape(tracks(1,2,:),[],1),...
            reshape(disptracks(f,1,:),[],1),...
            reshape(disptracks(f,2,:),[],1),...
            0,'-m');
    %plot(VxxPTV(isnan(Vqxn(:,:,f))),VyyPTV(isnan(Vqxn(:,:,f))),'ow');
    
    %plot displacement field
    %quiver(Vxx,Vyy,DSCALE*Vqx(:,:,f),DSCALE*Vqy(:,:,f),0,'-y');
    
    %plot calculated force
    %quiver(Vxx,Vyy,FSCALE*SF(:,:,1,f),FSCALE*SF(:,:,2,f),0,'-y');
    
    
    %Time Stamp
    if SHOW_TIMESTAMP
        str = sprintf('Time: %04.01f min',(Time(f)-Time(1))/60);
        text(xloc,yloc,str,'Color','w','VerticalAlignment','top','HorizontalAlignment','Left','FontSize',TIME_FONT_SIZE);
    end
    
    %Plot ScaleBar
    text(mean(SB_X),mean(SB_Y)+2*SB_WIDTH,sprintf('%d µm',SB_LENGTH),'color','w','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',SB_FONT_SIZE);
    plot(SB_X,SB_Y,'-w','LineWidth',SB_WIDTH);
    
    %save frame to animation
    Anim(f) = getframe(hax);
end
try
close(hfig);
catch
end

%% Save Movie
writerObj = VideoWriter(fullfile(mov_path,mov_file),'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);
writeVideo(writerObj,Anim);
close(writerObj); 