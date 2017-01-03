[f,d] = uigetfile();

TFMForceViewer(fullfile(d,f),'CellImageCLim',[650,1600],'SMAGLim',[0,3000],'PlotDisplacements',true);