num_points = 200;
num_entities = 20;
%num_entities = 201;
entity_width = 2/num_entities;
[X,Y,Z] = meshgrid(linspace(-1,1,num_points)); % mesh for calculation of exact solution
[Xc,Yc,Zc] = meshgrid(linspace(-1.+entity_width/2,1-entity_width/2,num_entities)); % mesh consisting of entity centers for projection
%[X,Y,Z] = meshgrid(linspace(-1,1,num_entities));
%psi0fun = @(X,Y,Z) max(1/(8*pi*0.03^2)*exp(-(X.^2+Y.^2+Z.^2)./2/0.03^2),1e-4/4/pi);
psi0fun = @(X,Y,Z) max(1/(0.03^3*pi^3)*exp(-1*(X.^2+Y.^2+Z.^2)./(pi*0.03^2)),1e-4/4/pi);
y = Pointsource(psi0fun,X,Y,Z,1);
Xcr = reshape(Xc, [num_entities^3 1]);
Ycr = reshape(Yc, [num_entities^3 1]);
Zcr = reshape(Zc, [num_entities^3 1]);
Xr = reshape(X, [num_points^3 1]);
Yr = reshape(Y, [num_points^3 1]);
Zr = reshape(Z, [num_points^3 1]);
yr = reshape(y, [num_points^3 1]);
ycr = zeros(num_entities^3, 1);

for ii = 1:num_entities^3
    indices = find(abs(Xr-Xcr(ii)) < entity_width/2+1e-10);
    indices = find(abs(Yr(indices)-Ycr(ii)) < entity_width/2+1e-10);
    indices = find(abs(Zr(indices)-Zcr(ii)) < entity_width/2+1e-10);
    [m,n] = size(indices);
    ycr(ii) = sum(yr(indices))/m;
    ii
end 
% renormalize
ycrnorm = ycr / (sum(ycr*((2./num_entities)^3)));
yrnorm = yr / (sum(yr*((2./num_points)^3)));

fileID = fopen('values_matlab.txt', 'w');
fprintf(fileID,'[%15.15f %15.15f %15.15f] \t %15.15f\n', [Xcr'; Ycr'; Zcr'; ycrnorm']);
fclose(fileID);

yindices = find(abs(Yr)<entity_width/19.);
zindices = find(abs(Zr)<entity_width/19.);
ycindices = find(abs(Ycr-entity_width/2.)<1e-10);
zcindices = find(abs(Zcr-entity_width/2.)<1e-10);
%yindices = find(abs(Yr)<1e-10);
%zindices = find(abs(Zr)<1e-10);
indices = intersect(yindices, zindices);
cindices = intersect(ycindices, zcindices);
plot(Xr(indices), yrnorm(indices), Xcr(cindices), ycrnorm(cindices))