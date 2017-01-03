function [Sx,Sy,Skx,Sky,kxx,kyy,M] = disp2stressFTTC(Ux,Uy,dx,dy,E,v)
% Calculates the Stress-Field for a given displacement using the Fourier
% Space solution to the inverse Boussinesq problem (shear stress on an
% infinite half-space).
%
% Inputs:
%   Ux,Uy:  The strain-fields imposed on the top of the elastic surface
%           These scalar-fields should be regularly space, with distance
%           between elements [dy,dx]
%
%   dx,dy:  The distance between points in the strain field
%           Should be in the units of Ux,Uy
%
%   E:     Young's modulous of the material
%   v:      Poisson's Ratio for the material
%           
% Output:
%   S:  2-D stress vector field
%           Where (nRows,nCols) = size(Ux)
%           (nRows,nCols,2) = size(S)
%           S(:,:,1) is the x-direction, S(:,:,2) is the y direction
%   Skx,Sky: Fourier transform of Stress field generated during
%            calculation. Useful for calculating stress moments
%   kxx,kyy: wave-vectors overwhich Skx,Sky are evaluated
% ************************************************************************
% This calculation is based on the FTTC methods presented in"
% Butler, et al. "Traction fields, moments and strain energy that cells
% exert on their surroundings." American Journal of Physiology: Cell
% Physiology, (2002) Vol 282, Iss. 3 pp 595-605

%% Changes Log: ===========================================================
%   2014-11-21: Initial creation by Daniel T. Kovari
%               Georgia Inst. of Tech, School of Physics
%               www.curtisresearch.gatech.edu
%==========================================================================

%% Check Syntax
if any(size(Ux)~=size(Uy))
    error('Ux and Uy must be the same size');
end

[Ny,Nx] = size(Ux);

%Calculate rescale factor so that matrix op work near integers
RS = sqrt(dx^2+dy^2);
dx = dx/RS;
dy = dy/RS;
Ux = Ux/RS;
Uy=Uy/RS;

Wsx = 2*pi/dx;
Wsy = 2*pi/dy;
kx = (0:(Nx-1))*(Wsx/Nx);
kx(kx>=Wsx/2) = kx(kx>=Wsx/2)-Wsx;
ky = (0:(Ny-1))*(Wsy/Ny);
ky(ky>=Wsy/2) = ky(ky>=Wsy/2)-Wsy;

[kxx,kyy] = meshgrid(kx,ky);

%Calculate inverse to get force
%Calc FT(disp)
Ukx = fft2(Ux);
Uky = fft2(Uy);
%calculate G^-1
iGK = zeros(Ny,Nx,4); %holder for G^-1

for k = 1:Nx
    for l = 1:Ny
        tmp = InvBoussinesqGk([kxx(l,k),kyy(l,k)],1,v);%put E=1 here to help with scaling
        %handle problem with indistinguishable nyquist frequency components
        if kxx(l,k)==pi/dx||kyy(l,k)==pi/dy
            %'here'
            tmp(1,2)=0;
            tmp(2,1)=0;
        end
        iGK(l,k,1) = tmp(1,1);
        iGK(l,k,2) = tmp(1,2);
        iGK(l,k,3) = tmp(2,1);
        iGK(l,k,4) = tmp(2,2);
    end
end

%calc FT(F)
Skx = iGK(:,:,1).*Ukx+iGK(:,:,2).*Uky;
Sky = iGK(:,:,3).*Ukx+iGK(:,:,4).*Uky;

%figure(99);
%imagesc(abs(Skx));
%Skx(kxx==0&kyy==0)

%inverse FT to get stress-field
Sx = ifft2(Skx,'symmetric')*E;%rescale result here 
Sy = ifft2(Sky,'symmetric')*E;

if nargout>2
    Skx=Skx*E;
    Sky=Sky*E;
end

if nargout>4
    kxx = kxx/RS;
    kyy=kyy/RS;
end

if nargout>6 %compute strain moment
    x = (0:Nx-1)*dx*RS;
    %dx*RS
    %max(Sx(:))
    y = (0:Ny-1)*dy*RS;
    [xx,yy] = meshgrid(x,y);
    M = zeros(2);
    M(1,1) = sum(sum(xx.*Sx));
    M(1,2) = 1/2*sum(sum(xx.*Sy +yy.*Sx));
    M(2,1) = 1/2*sum(sum(xx.*Sy +yy.*Sx));
    M(2,2) = sum(sum(yy.*Sy));
    %M
end



% function tiG = TikhonovInverseG(K,E,v,lambda)
% k = sqrt(K(1)^2+K(2)^2);
% G = 2*(1+v)/E*(1/k^3)*[k^2-v*K(1)^2, -v*K(1)*K(2);...
%                         -v*K(1)*K(2), k^2-v*K(2)^2];
% % M = G'*G+lambda*eye(2);
% 
% tiG = (G'*G+lambda*eye(2))^-1*G';
% % d = M(1,1)*M(2,2)-M(1,2)*M(2,1)
% % C = [M(2,2),-M(1,2);-M(2,1),M(1,1)]
% % tiG = (1/d)*C*G';

function invGk = InvBoussinesqGk(K,E,v)

k = sqrt(K(1)^2+K(2)^2);

if k==0
    %k=eps;
    invGk = [0,0;0,0];
    return;
end

invGk =E/(2*(1-v^2)*k)*...
    [k^2-v*K(2)^2, v*K(1)*K(2);...
     v*K(1)*K(2), k^2-v*K(1)^2];

