clc;
IM1 = rand(1000,1000,5);
I = imdilate(eye(1000,1000),ones(5));
I = repmat(I,1,1,5);

[hFig,hAx,hCB] = overlay_animfig(IM1,randn(1000,1000,5),'CLim','scaled', 'colormap',gray(256));

ylabel(hCB,'test');