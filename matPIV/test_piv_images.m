clc;

im2 = bpass((imread('sample002_1_1.tif')),1,7);
im1 = bpass((imread('sample002sds_1_1.tif')),1,7);

[H,W] = size(im1);

tic
[vx,vy,x,y,vx_nan,vy_nan] = piv_rec(mat2gray(im1),mat2gray(im2),...
    'startW',64,...
    'endW',32,...
    'startH',64,...
    'endH',32);
toc

figure(1);
clf;
image('cdata',repmat(reshape([1,0,0],1,1,3),H,W,1),'alphadata',im1,'alphadatamapping','scaled');
axis image;

set(gca,'color','k');
hold on;

image('cdata',repmat(reshape([0,1,1],1,1,3),H,W,1),'alphadata',im2,'alphadatamapping','scaled');

[xx,yy] = meshgrid(x,y);

vx = vx - nanmean(vx(:));
vy = vy -  nanmean(vy(:));

[vx,vy,qxn,qyn] = residfilt(vx,vy);

quiver(xx,yy,vx,vy,0);

xn = xx(isnan(qxn));
yn = yy(isnan(qyn));
plot(xn(:),yn(:),'sw');

