im = mat2gray(im2double(imread('t001xy1c1.tif')));


sz = 9;
B = bpass(im,1,sz);

figure(1); clf;
subplot(1,2,1);
imagesc(im); axis image; colormap jet; set(gca,'CLim',[0,0.2]);hold on;
subplot(1,2,2);
imagesc(B); axis image; colormap jet; set(gca,'CLim',[0,0.2]);hold on;

TH = 0.055;

clc;

disp('old pkfnd method');
tic
[PKS,IND] = pkfnd(B,0.02,sz);
toc
CNT=pkfnd_subpixel(im,PKS,'centroid');

%CNT = cntrd(B,PKS,sz);
figure(1);
subplot(1,2,2);
plot(PKS(:,1),PKS(:,2),'sg');
plot(CNT(:,1),CNT(:,2),'+g');
disp(['Number of CNT: ',num2str(length(CNT))]);

disp('new pkfnd method');
%tic

%[PKS2,IND2] = pkfnd2(B,sz,mean(B(:)));
%toc
%CNT2 = cntrd(B,PKS2,sz);
%CNT2=pkfnd_subpixel(im,PKS2,'centroid');

[XY,SNR] = imfindpeak(im)

figure(1);
subplot(1,2,2);
%plot(PKS2(:,1),PKS2(:,2),'om');
%plot(CNT2(:,1),CNT2(:,2),'+m');
plot(XY(1),XY(2),'sy','markersize',24);
%disp(['Number of CNT2: ',num2str(length(CNT2))]);

