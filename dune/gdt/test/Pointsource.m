function y = Pointsource(psi0fun,X,Y,Z,t,rhoD,fD)
% Example: 
%[X,Y,Z] = meshgrid(linspace(-1,1));
%psi0fun = @(X,Y,Z) max(1/sqrt(8*pi*0.03^2)*exp(-(X.^2+Y.^2+Z.^2)./2/0.03^2),1e-4/4/pi);
%y = Pointsource(psi0fun,X,Y,Z,1);
%contourslice(X,Y,Z,y,-1:0.25:1,-1:0.25:1,-1:0.25:1);
R = sqrt(X.^2+Y.^2+Z.^2);

%check 1D mass
rhoD = rhoD(1:end-2);
fD = fD(1:end-2);
dx = diff(rhoD);
Mass1D = sum((fD(2:end)+fD(1:end-1))/2.*dx)
%tD = fD.*rhoD.^2;
%Mass3D = sum((tD(2:end)+tD(1:end-1))/2.*dx)
dx3 = diff(rhoD.^3);
Mass3D = sum((fD(2:end)+fD(1:end-1))/2.*dx3)*4/3*pi

funD = @(r) interp1(rhoD,fD,r,'spline',0);

rho = funD(R);
%plot(R(:), rho(:))

psi0 = psi0fun(X,Y,Z); %Smooth initial condition, assume size(X) = [ny,nx,nz];


% calculate integrals of psi0 and rho
initialmass = sum(psi0fun(X(:),Y(:),Z(:))/numel(X))*(max(X(:))-min(X(:)))*(max(Y(:))-min(Y(:)))*(max(Z(:))-min(Z(:)))
massdirac = sum(rho(:)./numel(X))*(max(X(:))-min(X(:)))*(max(Y(:))-min(Y(:)))*(max(Z(:))-min(Z(:)))
%initialmass =
%integral3(psi0fun,min(X(:)),max(X(:)),min(Y(:)),max(Y(:)),min(Z(:)),max(Z(:)));
%y = convn(rho,psi0,'same');
y = convolution3D_FFTdomain(rho,psi0);
y = y./numel(y)*(max(X(:))-min(X(:)))*(max(Y(:))-min(Y(:)))*(max(Z(:))-min(Z(:))); %map to quadrature
finalmass = sum(y(:)/numel(X))*(max(X(:))-min(X(:)))*(max(Y(:))-min(Y(:)))*(max(Z(:))-min(Z(:)))
finalmassexpected = massdirac*initialmass % it should be int(f*g) = int(f)*int(g), but appears to be 0.5*int(g) instead...
factorfinalmass = finalmassexpected/finalmass


end

