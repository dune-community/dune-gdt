function y = Pointsource(psi0fun,X,Y,Z,t)
% Example: 
%[X,Y,Z] = meshgrid(linspace(-1,1));
%psi0fun = @(X,Y,Z) max(1/sqrt(8*pi*0.03^2)*exp(-(X.^2+Y.^2+Z.^2)./2/0.03^2),1e-4/4/pi);
%y = Pointsource(psi0fun,X,Y,Z,1);
%contourslice(X,Y,Z,y,-1:0.25:1,-1:0.25:1,-1:0.25:1);
R = sqrt(X.^2+Y.^2+Z.^2);

[rhoD,fD] = point_source_solution(t); %1D solution
%check 1D mass

rhoD = rhoD(1:end-2);
fD = fD(1:end-2);
dx = diff(rhoD);
Mass1D = sum((fD(2:end)+fD(1:end-1))/2.*dx)
tD = fD.*rhoD.^2;
Mass3D = sum((tD(2:end)+tD(1:end-1))/2.*dx)

funD = @(r) interp1(rhoD,fD,r,'linear',0);

rho = funD(R);

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

function [rho,f] = point_source_solution(t)
% Evaluates the true solution of the line source test case,
% using Ganapol's quadrature formulas.
rho = (1-linspace(1,0,1000).^2)*t;
f = phi_pt(rho,t)*pi;
end

function f = phi_pt(r,t)
r = max(r,1e-10); % move zero value into small positive regime
eta = r/t;
ind = find(eta<1);
g = r*0;
% compute quad
for k = ind
    g(k) = integral(@(u) sec(u/2).^2.*real(...
        (eta(k)+1i*tan(u/2)).*...
        xi(eta(k),u).^3.*exp(t/2*(1-eta(k).^2).*xi(eta(k),u)) ...
        ),0,pi);
end
% final result
f = 1/(2*pi)*exp(-t)./(4*pi*r.*t.^2)*(t/2)^2.*(1-eta.^2).*g;
f(ind) = f(ind)+phi_pt0(r(ind),t)+phi_pt1(r(ind),t);
end

function f = xi(eta,u)
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
f = (log(q)+1i*u)./(eta+1i*tan(u/2));
end

function y = phi_pt0(r,t)
diractol = 1e-8;
eta = r/t;
diracev = 1/sqrt(2*pi*diractol)*exp(-(1-eta).^2/2/diractol);
y = exp(-t)/4/pi./r./t.^2.*diracev;
end

function y = phi_pt1(r,t)
eta = r/t;
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
y = exp(-t)./(4*pi*r*t^2)*t.*log(q);
end