IM1 = rand(100,100,5);
I = imdilate(eye(100,100),ones(5));
I(~I) = NaN;
I = repmat(I,1,1,5);

[hFig,hAx,hImages,hCB] = composit_animfig(IM1,I,'CLim','scaled', 'ShowColorbar',true,'ColorbarWidth',20,'colormap',jet(256));

set(hFig,'units','inches','position',[0,0,5,5]);

ylabel(hCB,'test');

