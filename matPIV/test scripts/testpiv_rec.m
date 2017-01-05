im1 = mat2gray(imread('t001xy1c1.tif'));
im2 = mat2gray(imread('t124xy1c1.tif'));
im1 = bpass(im1,1,11);
im2 = bpass(im2,1,11);

% clc;
% rng(5000);
% clear all;
% %construct image
% W=512;
% H=512;
% dXm = 0;
% dYm = 0;
% 
% driftX = -5;
% driftY = 5;
% 
% sig2 = 3;
% nP = 2000;
% 
% XY = [rand(nP,1)*(W-1)+1,rand(nP,1)*(H-1)+1];
% dXY = [(2*rand(nP,1)-1)*dXm+driftX,(2*rand(nP,1)-1)*dYm+driftY];
% 
% mean(dXY,1)
% 
% im1=zeros(H,W);
% im2=zeros(H,W);
% for p=1:nP
%     im1=im1+imgaussian([sig2,sig2],XY(p,:),[H,W]);
%     im2=im2+imgaussian([sig2,sig2],XY(p,:)+dXY(p,:),[H,W]);
% end

figure(1); clf;
subplot(1,2,1);
imagesc(im1);colormap jet; axis image;
subplot(1,2,2);
imagesc(im2);colormap jet; axis image;
drawnow;

[vx,vy,x,y,vx_nan,vy_nan] = piv_rec(im1,im2,'startH',32,'startW',32,'endH',16,'endW',16);

%eliminate drift
vx2 = vx-nanmedian(vx_nan(:));
vy2 = vy-nanmedian(vy_nan(:));

%fix precision to 0.2
vx2 = 0.2*fix(vx2/0.2);
vy2 = 0.2*fix(vy2/0.2);

[xx,yy] = meshgrid(x,y);

%clamp edges to zero
vx2(xx<32|xx>size(im1,2)-32) = 0;
vy2(xx<32|xx>size(im1,2)-32) = 0;
vx2(yy<32|yy>size(im1,1)-32) = 0;
vy2(yy<32|yy>size(im1,1)-32) = 0;

figure(1);
subplot(1,2,1);
hold on;
quiver(xx,yy,10*vx2,10*vy2,0,'-g');
plot(xx(isnan(vx_nan)),yy(isnan(vx_nan)),'oc');
title('vectors scaled by 10x, o specifies poorly defined vectors');