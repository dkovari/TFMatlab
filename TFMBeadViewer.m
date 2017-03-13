function TFMBeadViewer(TFMdata,varargin)

p = inputParser;
p.CaseSensitive = false;
addParameter(p,'MoviePath',[]);
addParameter(p,'CloseAfterSave',true);
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
fields = {'cnt','tracks','imstack','disptracks','Time','Ostack'};
if any( ~isfield(TFMdata,fields))
    erfld = fields(~isfield(TFMdata,fields));
    error('TFMdata Missing Field: %s',erfld{:});
end

[hFig,hAx,hCB] = overlay_animfig(TFMdata.Ostack,TFMdata.Bstack,...
    'OverlayColor',[0,1,1],...
    'CLim','average',...
    'colormap',gray(256),...
    'frameupdate_fn',@FrameChange);

ylabel(hCB,'Bead Intensity');

hold(hAx,'on');
hCnt = plot(hAx,TFMdata.cnt{2}(:,1), TFMdata.cnt{2}(:,2),'+k','markersize',10);
hQvr = quiver(reshape(TFMdata.tracks(1,1,:),[],1),...
            reshape(TFMdata.tracks(1,2,:),[],1),...
            reshape(TFMdata.disptracks(1,1,:),[],1),...
            reshape(TFMdata.disptracks(1,2,:),[],1),...
            0,'-m');

    function FrameChange(~,~,f)
        set(hCnt,'XData',TFMdata.cnt{f+1}(:,1),...
                 'YData',TFMdata.cnt{f+1}(:,2));
        set(hQvr,'UData',reshape(TFMdata.disptracks(f,1,:),[],1),...
                 'VData',reshape(TFMdata.disptracks(f,2,:),[],1));
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