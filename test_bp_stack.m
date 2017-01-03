bpass_lnoise = 1;   %approximate size of noise, in pixles (usually 1)
bpass_sz = 9;       %approximate particle size, in pixels. objects larger and smaller will be filtered with a spatial bandpass
pkfnd_sz = 9;       %approx. particle size, in pixels. size of the window used to find particle center of mass
pkfnd_th = 1;    %minimum peak intensity after bandpass operation
cnt_sz = 9;         %particle size, in pixels

nF = size(imstack,3);
imstack(1:50,:,:) = NaN;
Bstack = zeros(size(imstack));

for f=1:nF
    Bstack(:,:,f) = bpass(imstack(:,:,f),bpass_lnoise,bpass_sz);
end