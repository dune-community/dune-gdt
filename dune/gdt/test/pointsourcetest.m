num_points = 501;
num_entities = 30;
t_end = 0.5
entity_width = 2/num_entities;
[X,Y,Z] = meshgrid(linspace(-1,1,num_points)); % mesh for calculation of exact solution
[Xc,Yc,Zc] = meshgrid(linspace(-1.+entity_width/2,1-entity_width/2,num_entities)); % mesh consisting of entity centers for projection
%psi0fun = @(X,Y,Z) max(1/(8*pi*0.03^2)*exp(-(X.^2+Y.^2+Z.^2)./2/0.03^2),1e-4/4/pi);
psi0fun = @(X,Y,Z) max(1/(4*0.03^3*pi^4)*exp(-1*(X.^2+Y.^2+Z.^2)./(pi*0.03^2)),1e-4/4/pi);
y = Pointsource(psi0fun,X,Y,Z, t_end);
radius = entity_width/2+1e-10;
[Xcsorted, I] = sort(Xc(:));
Ycsorted = Yc(I);
Zcsorted = Zc(I);
ycsorted = zeros(num_entities^3, 1);
indices1 = find(abs(X(:)-Xcsorted(1)) < radius);
indices2 = find(abs(Y(indices1)-Ycsorted(1)) < radius);
lastX = Xcsorted(1);
lastY = Ycsorted(1);
 for ii = 1:num_entities^3
     if Xcsorted(ii) ~= lastX
       indices1 = find(abs(X(:)-Xcsorted(ii)) < radius);
     end
     if (Xcsorted(ii) ~= lastX) || (Ycsorted(ii) ~= lastY)
        indices2 = find(abs(Y(indices1)-Ycsorted(ii)) < radius);
            lastY = Ycsorted(ii);
     end
     indices3 = find(abs(Z(indices1(indices2))-Zcsorted(ii)) < radius);
     indices = indices1(indices2(indices3));
     [m,n] = size(indices);
     ycsorted(ii) = sum(y(indices))/m;
     lastX = Xcsorted(ii);
     if mod(ii,1000) == 0
         ii
     end
 end 

fileID = fopen(string('values_matlab_') + string(num_entities) + string('.txt'), 'w');
fprintf(fileID,'[%15.15f %15.15f %15.15f] \t %15.15f\n', [Xcsorted(:)'; Ycsorted(:)'; Zcsorted(:)'; ycsorted(:)']);
fclose(fileID);

ycindices = find(abs(Ycsorted-entity_width/2.)<1e-10);
zcindices = find(abs(Zcsorted-entity_width/2.)<1e-10);
yindices = find(abs(Y(:))<1e-10);
zindices = find(abs(Z(:))<1e-10);
indices = intersect(yindices, zindices);
cindices = intersect(ycindices, zcindices);
plot(X(indices), y(indices), Xcsorted(cindices), ycsorted(cindices))
print('pointsourceplot500','-dpng')
contourslice(X,Y,Z,y,-1:0.25:1,-1:0.25:1,-1:0.25:1)
