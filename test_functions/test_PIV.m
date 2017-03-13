%Test 
[H,W] = size(TFMdata.Iref);

figure(1);
clf;

image('cdata',repmat(reshape([1,0,0],1,1,3),H,W,1),'alphadata',TFMdata.Iref);
axis image;

set(gca,'color','k');
hold on;

image('cdata',repmat(reshape([0,1,1],1,1,3),H,W,1),'alphadata',TFMdata.Bstack(:,:,1));

[vx,vy,x,y,vx_nan,vy_nan] = piv_rec(TFMdata.Iref,TFMdata.Bstack(:,:,1),...
    'startW',100,...
    'endW',20,...
    'startH',100,...
    'endH',20);

[xx,yy] = meshgrid(x,y);
quiver(xx,yy,vx,vy,0);