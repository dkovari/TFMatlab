f = 1;
d = [];
FilePath = {};
while f~=0
    [f,d] = uigetfile([d,'*.mat']);
    if f~=0
        FilePath = [FilePath,fullfile(d,f)];
    end
end

%%
for FP = FilePath
    [d,n,e] = fileparts(FP{1});
    TFMForceViewer(FP{1},'CellImageCLim',[650,1600],...
        'SMAGLim',[0,3000],...
        'PlotStrain',true,...
        'MoviePath',fullfile(d,[n,'.mp4']),...
        'FigureSize',[1080,720]);
end