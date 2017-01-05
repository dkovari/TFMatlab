s=32;
e=55.25;
dx = e-s;
im1=double(BWcircle(4,[s,16],[32,64]));
im2=double(BWcircle(4,[e,16],[32,64]));

[dXY,SNR,C]=imxcorr(im1,im2);
disp(dXY)
disp(SNR)
figure(3);clf;
subplot(2,2,1);
imagesc(im1);axis image; colormap jet;
subplot(2,2,2);
imagesc(im2);axis image; colormap jet;
title(sprintf('dx=%0.2f',dx));
subplot(2,2,[3,4]);
imagesc(C); axis image; colormap jet;
hold on;
quiver(64,32,dXY(1),dXY(2),0,'+-k');
plot(64+dXY(1),32+dXY(2),'xk');
title(sprintf('measured dx=%0.2f',dXY(1)));
