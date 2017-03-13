Dirs = {};
FilePath = {};
d = 1;
LastDir = '';
while d~=0
    d = uigetdir(LastDir);
    if d~=0
        Dirs=[Dirs,d];
    end
    LastDir = d;
end

%%
for j = 1:numel(Dirs)
    fprintf('Animating %s (%d/%d):\n',Dirs{j},j,numel(Dirs));
    files = dir(fullfile(Dirs{j},'*.mat'));
    for f=1:numel(files)
        fprintf('\t %d/%d\n',f,numel(files));
        [~,n,e] = fileparts(files(f).name);
        TFMForceViewer(fullfile(Dirs{j},files(f).name),...
            'CellImageCLim',[650,1600],...
            'PlotStrain',true,...
            'MapData','StrainEnergyDensity',...
            'MapLim',[0,6e-4],...
            'MoviePath',fullfile(Dirs{j},[n,'.mp4']),...
            'MovieQuality',75,...
            'FigureSize',[1280,1010]);
    end
end