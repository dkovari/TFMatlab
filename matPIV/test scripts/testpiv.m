% im1 = mat2gray(imread('t001xy1c1.tif'));
% im2 = mat2gray(imread('t124xy1c1.tif'));
% im1 = bpass(im1,1,11);
% im2 = bpass(im2,1,11);
clc;
close all;
rng(5000);
clear all;
%construct image
W=128;
H=128;
dXm = 0;
dYm = 0;

driftX = 10;
driftY = 0;

sig2 = 3;
nP = 200;

XY = [rand(nP,1)*(W-1)+1,rand(nP,1)*(H-1)+1];
dXY = [(2*rand(nP,1)-1)*dXm+driftX,(2*rand(nP,1)-1)*dYm+driftY];

mean(dXY,1)

im1=zeros(H,W);
im2=zeros(H,W);
for p=1:nP
    im1=im1+imgaussian([sig2,sig2],XY(p,:),[H,W]);
    im2=im2+imgaussian([sig2,sig2],XY(p,:)+dXY(p,:),[H,W]);
end

figure(1); clf;
subplot(1,2,1);
imagesc(im1);colormap jet; axis image; set(gca,'CLim',[0,0.03]);
subplot(1,2,2);
imagesc(im2);colormap jet; axis image; set(gca,'CLim',[0,0.03]);

%mean drift

C = xcorr2(im2,im1);
figure(2); clf;
imagesc(C);colormap jet; axis image;

pks = pkfnd2(C,5);
IND = sub2ind(size(C),pks(:,2),pks(:,1));
[mxpkv,mxpk] = max(C(IND));

%find second highest peak
if size(pks,1)>2
    pks2 = pks;
    pks2(mxpk,:) = [];
    IND = sub2ind(size(C),pks2(:,2),pks2(:,1));
    [mxpkv2,mxpk2] = max(C(IND));
elseif size(pks,1)==2
    if mxpk == 1
        mxpk2 = 2;
    else
        mxpk2=1;
    end
    mxpkv2 = C(pks2(mxpk2,1),pks2(mxpk2,1));
else
    mxpk2=NaN;
    mxpk2v=0;
end
%calc SNR (need to look this up)
SNR = mxpkv/mxpkv2

XY = pkfnd_subpixel(C,pks(mxpk,:),'gaussian')

figure(2); hold on;
plot(XY(1),XY(2),'+w');
plot(pks(:,1),pks(:,2),'sw');

meandrift = XY-[W,H]

%begin recursive PIV

vx = repmat(meandrift(1),3);
vy = repmat(meandrift(2),3);
x=[0.5,64.5,128.5];
y=[0.5,64.5,128.5];

[vx2,vy2,x2,y2]=piv(im1,im2,64,64,x,y,vx,vy,'interp');


[xx,yy] = meshgrid(x2,y2);

figure(1);
subplot(1,2,1);
hold on;
quiver(xx,yy,vx2,vy2,0,'-g');
drawnow;

return;

[vx2,vy2,x2,y2]=piv(im1,im2,32,32,x2,y2,vx2,vy2,'interp');


[xx,yy] = meshgrid(x2,y2);

figure(1);
subplot(1,2,1);
hold on;
quiver(xx,yy,vx2,vy2,0,'-m');
drawnow;
[vx2,vy2,x2,y2]=piv(im1,im2,16,16,x2,y2,vx2,vy2,'interp');


[xx,yy] = meshgrid(x2,y2);

figure(1);
subplot(1,2,1);
hold on;
quiver(xx,yy,vx2,vy2,0,'-c');
drawnow;

disp('mean piv')
disp([mean(vx2(:)),mean(vy2(:))]);




