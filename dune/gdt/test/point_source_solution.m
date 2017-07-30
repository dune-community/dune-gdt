function [rho,f] = point_source_solution(t)
% Evaluates the true solution of the line source test case,
% using Ganapol's quadrature formulas.
rho = (1-linspace(1,0,100000).^2);
f = phi_pt(rho,t);
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
%f = 1/(2*pi)*exp(-t)./(4*pi*r.*t.^2)*(t/2)^2.*(1-eta.^2).*g;
%f(ind) = f(ind)+phi_pt0(r(ind),t)+phi_pt1(r(ind),t);
f = phi_pt0(r,t);
f(ind) = f(ind) + 1/(2*pi)*exp(-t)./(4*pi*r(ind).*t.^2)*(t/2)^2.*(1-eta(ind).^2).*g(ind) + phi_pt1(r(ind),t);
end

function f = xi(eta,u)
q = (1+eta)./(1-eta);
q = min(max(q,-1e80),1e80); % remove singularity at eta = 1
f = (log(q)+1i*u)./(eta+1i*tan(u/2));
end

function y = phi_pt0(r,t)
diractol = 1e-6;
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