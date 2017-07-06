function y = Pointsource(psi0fun,X,Y,Z,t)
% Example: 
%[X,Y,Z] = meshgrid(linspace(-1,1));
%psi0fun = @(X,Y,Z) max(1/sqrt(8*pi*0.03^2)*exp(-(X.^2+Y.^2+Z.^2)./2/0.03^2),1e-4/4/pi);
%y = Pointsource(psi0fun,X,Y,Z,1);
%contourslice(X,Y,Z,y,-1:0.25:1,-1:0.25:1,-1:0.25:1);
R = sqrt(X.^2+Y.^2+Z.^2);

[rhoD,fD] = line_source_solution(t); %1D solution
rhoD = rhoD(1:end-2);
fD = fD(1:end-2);
funD = @(r) interp1(rhoD,fD,r,'linear',0);

rho = funD(R);

psi0 = psi0fun(X,Y,Z); %Smooth initial condition, assume size(X) = [ny,nx,nz];

%y = convn(rho,psi0,'same');
y = convolution3D_FFTdomain(rho,psi0);
y = y/numel(y); %map to quadrature




end

function [rho,f] = line_source_solution(t)
% Evaluates the true solution of the line source test case,
% using Ganapol's quadrature formulas.
rho = [(1-linspace(1,0,200).^2)*t 1];
f = phi_l(rho,t)*pi;
end

function f = phi_l(rho,t)
eta = rho/t;
ind = find(eta<1);
f = rho*0; phi_l0 = f;
for k = ind
    f(k) = quad(@(w) phi_pt(t*sqrt(eta(k).^2+w.^2),t),...
    0,sqrt(1-eta(k).^2),1e-5);
end
phi_l0(ind) = exp(-t)/(2*pi*t^2)./sqrt(1-eta(ind).^2);
f = phi_l0+(2*t)*f;
end

function f = phi_pt(r,t)
r = max(r,1e-10); % move zero value into small positive regime
eta = r/t;
ind = find(eta<1);
g = r*0;
% compute quad
for k = ind
    g(k) = quad(@(u) sec(u/2).^2.*real(...
        (eta(k)+1i*tan(u/2)).*...
        xi(eta(k),u).^3.*exp(t/2*(1-eta(k).^2).*xi(eta(k),u)) ...
        ),0,pi,1e-2);
end
% final result
f = 1/(2*pi)*exp(-t)./(4*pi*r.*t.^2)*(t/2)^2.*(1-eta.^2).*g;
f(ind) = f(ind)+phi_pt1(r(ind),t);
end

function f = xi(eta,u)
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
f = (log(q)+1i*u)./(eta+1i*tan(u/2));
end

function y = phi_pt1(r,t)
eta = r/t;
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
y = exp(-t)./(4*pi*r*t^2)*t.*log(q);
end