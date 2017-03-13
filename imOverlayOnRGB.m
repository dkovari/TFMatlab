function O = imOverlayOnRGB(I_baseRGB,I_over,RGB_over,clim_in,clim_out)

%recast clim to be the same type as I_over

clim_in = cast(clim_in,class(I_over));
clim_out = cast(clim_out,class(I_over));
%shift and scale the overlay values to alpha range
I_over = mat2gray(I_over,clim_in).*(clim_out(2)-clim_out(1))+clim_out(1);

src = cat(3,RGB_over(1)*I_over,RGB_over(2)*I_over,RGB_over(3)*I_over);

dst = cat(3, I_baseRGB(:,:,1).*(1-I_over(:,:)) ,...
             I_baseRGB(:,:,2).*(1-I_over(:,:)) ,...
             I_baseRGB(:,:,3).*(1-I_over(:,:)) );

O = imadd(src,dst);