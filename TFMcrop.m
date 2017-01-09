function TFMdata2 = TFMcrop(TFMdata,varargin)
%Interactive Tool for cropping TFM data
%
% Input:
%   Accepts same arguments as TFMForceViewer

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
    [Dir,Name,~] = fileparts(TFMdata);
    hDlg = msgbox({'Loading Data','Please wait, this will take a while.'},'Loading...');
    TFMdata = load(TFMdata);
    try
        delete(hDlg);
    catch
    end
else
    Dir = [];
    Name = '';
end
%% Validate TFMdata
fields = {'Time','Ostack','Vxx','Vyy','Vqx','Vqy','SF','dx','dy','E','v','PX_SCALE','SMAG'};
if any( ~isfield(TFMdata,fields))
    erfld = fields(~isfield(TFMdata,fields));
    error('TFMdata Missing Field: %s\n',erfld{:});
end
%TFM viewer
[hFig,hAx,hFig_hist,hAx_hist] = TFMForceViewer(TFMdata,varargin{:});
%% Prompt user for Polygon
title(hAx,'Click to draw a selection polygon. Close when done');
hPoly = impoly(hAx);
PolyPos = [];
origFigDel = get(hFig,'DeleteFcn');
    function FigDel(h,e)
        PolyPos = hPoly.getPosition;
        try
            origFigDel(h,e);
            delete(h);
        catch
        end
    end
set(hFig,'DeleteFcn',@FigDel);
waitfor(hFig)

%% Prepare TFMdata2
TFMdata2 = TFMdata;

%% Convert Poly to mask, crop to limits
X_CROP = [floor(min(PolyPos(:,1))),ceil(max(PolyPos(:,1)))];
Y_CROP = [floor(min(PolyPos(:,2))),ceil(max(PolyPos(:,2)))];

PolyMask = poly2mask(PolyPos(:,1),PolyPos(:,2),size(TFMdata.Ostack,1),size(TFMdata.Ostack,2));
PolyMaskV = imresize(PolyMask,size(TFMdata.Vxx),'nearest');

%% Apply mask to displacement field
TFMdata2.Vqx = bsxfun(@times,TFMdata.Vqx,PolyMaskV);
TFMdata2.Vqy = bsxfun(@times,TFMdata.Vqy,PolyMaskV);

%% Store mask data in file
TFMdata2.PolyPos = PolyPos;

%% Prompt to recalc stress
% answer = questdlg('Do you want to recalculate the stress field?','Recalculate?','yes','no','yes');
% if strcmp(answer,'yes')
    %% Calculate Stress-Field

    % Calculate stress using FTTC
    hWait = waitbar(0,'Calculating Traction Force using FTTC');

    %Force parameters
    PX_SCALE = TFMdata2.PX_SCALE;

    dx = TFMdata2.dx;
    dy = TFMdata2.dy;
    [Ny,Nx,~] = size(TFMdata2.Vqx);
    nF = size(TFMdata2.Vqx,3);
    StressField = zeros(Ny,Nx,2,nF);
    for f=1:nF
        [StressField(:,:,1,f),StressField(:,:,2,f)] = disp2stressFTTC(PX_SCALE*TFMdata2.Vqx(:,:,f),PX_SCALE*TFMdata2.Vqy(:,:,f),dx,dy,TFMdata2.E,TFMdata2.v);
        waitbar(f/nF,hWait);
    end
    try
    close(hWait);
    catch
    end

    %% Calculate Stress mag
    Stress_mag = zeros(size(TFMdata2.Vqx));
    nF = size(TFMdata2.Vqx,3);
    for f=1:nF
        Stress_mag(:,:,f) = sqrt(real(StressField(:,:,1,f)).^2+real(StressField(:,:,2,f)).^2);
    end

    %% Calculate Strain-energy density
    StrainEnergyDensity = zeros(size(TFMdata2.Vqx));
    StrainEnergy = zeros(nF,1);
    for f=1:nF
        StrainEnergyDensity(:,:,f) = PX_SCALE*0.5*real(StressField(:,:,1,f).*TFMdata2.Vqx(:,:,f)+StressField(:,:,2,f).*TFMdata2.Vqy(:,:,f));
        StrainEnergy(f) = sum(sum(StrainEnergyDensity(:,:,f)))*dx*dy;
    end
    
    %% Set Values
    TFMdata2.SF = StressField;
    TFMdata2.SMAG = Stress_mag;
    TFMdata2.SED = StrainEnergyDensity;
    TFMdata2.StrainEnergy = StrainEnergy;
% else
%     TFMdata2.Vqx = bsxfun(@times,TFMdata2.SF,PolyMaskV);
%     TFMdata2.SMAG = bsxfun(@times,TFMdata2.SF,PolyMaskV);
%     %% Calculate Strain-energy density
%     StrainEnergyDensity = zeros(size(TFMdata2.Vqx));
%     nF = size(TFMdata2.Vqx,3);
%     StrainEnergy = zeros(nF,1);
%     for f=1:nF
%         StrainEnergyDensity(:,:,f) = TFMdata2.PX_SCALE*0.5*real(TFMdata2.SF(:,:,1,f).*TFMdata2.Vqx(:,:,f)+TFMdata2.SF(:,:,2,f).*TFMdata2.Vqy(:,:,f));
%         StrainEnergy(f) = sum(sum(StrainEnergyDensity(:,:,f)))*TFMdata2.dx*TFMdata2.dy;
%     end
%     TFMdata2.SED = StrainEnergyDensity;
%     TFMdata2.StrainEnergy = StrainEnergy;
% end

%% Plot SE
hSE = figure();
plot((TFMdata.Time-TFMdata.Time(1))/60,StrainEnergy,'-r');

%moving average
M = movmean(StrainEnergy,5,'omitnan');
hold on;
hL = plot((TFMdata.Time-TFMdata.Time(1))/60,M,'--k');
title('Strain Energy');
ylabel('Strain Energy [J]');
xlabel('Time [min]');
legend(hL,'5-point moving average');

%% Prompt to save data
answer = questdlg('Do you want to save the data?','Save?','yes','no','yes');
if strcmp(answer,'yes')
    [FileName,PathName] = uiputfile(fullfile(Dir,[Name,'_crop.mat']),'Save file');
    if FileName~=0
        hDlg = msgbox({'Saving Data','Please wait, this will take a while.'},'Saving...');
        savefast(fullfile(PathName,FileName),...
            '-struct',TFMdata2);
        try
        close(hDlg);
        catch
        end
        
        [~,name,~] = fileparts(FileName);
        saveas(hSE,fullfile(PathName,[name,'_StrainEnergy.fig']))
        
    end
end

%% View Cropped Data
TFMForceViewer(TFMdata2,varargin{:});

%% Clear output variable if not requested
if nargout<1
    clear TFMdata2;
end

end
