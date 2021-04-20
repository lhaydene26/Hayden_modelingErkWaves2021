% model parameters
alpha_1 = 12.6;
beta_1 = 0.35;
gamma_1 = 15.4;
% alpha_2 = 0; % leave this uncommented to use the excited region
alpha_2 = 0.112; % leave this uncommented instead to use the source region
alpha_3 = 3.9;
gamma_2 = 11.8;
alpha_4 = .5; 
gamma_3 = 0.14;
gamma_e = 0.14;


% Define the input grid
lspace = 9127;
[x, y] = meshgrid(linspace(0, 0.4, lspace));
% Calculate the two surfaces
z1 = (gamma_2*x - alpha_2)/alpha_3;
z2 = y/alpha_4;
z3 = ((alpha_1*x.^2)./(beta_1^2 + x.^2))./(gamma_1*y + gamma_e + (alpha_1*x.^2)./(beta_1^2 + x.^2));

% Take the difference between each pairwise surface height and find the contour where those surfaces are zero.
zdiff1 = z1 - z2;
zdiff2 = z3 - z1;
zdiff3 = z3 - z2;
C1 = contours(x, y, zdiff1, [0 0]);
C2 = contours(x, y, zdiff2, [0 0]);
C3 = contours(x, y, zdiff3, [0 0]);
% Extract the x- and y-locations from the contour matrix C.
xL1 = C1(1, 2:end);
yL1 = C1(2, 2:end);
xL2 = C2(1, 2:end);
yL2 = C2(2, 2:end);
xL3 = C3(1, 2:end);
yL3 = C3(2, 2:end);
% Interpolate on the first surface to find z-locations for the intersection lines.
zL1 = interp2(x, y, z1, xL1, yL1);
zL2 = interp2(x, y, z1, xL2, yL2);
zL3 = interp2(x, y, z2, xL3, yL3);



% find intersections of each line
i = 1;
clear line_intersection_points
line_intersection_points = [];
line_intersection_points2 = [];
min_distance = 100;
min_xl = 0;
min_yl = 0;
min_zl = 0;
for j=1:numel(xL2)
    for k = 1:numel(xL3)
        dxL = abs(xL2(j) - xL3(k));
        dyL = abs(yL2(j) - yL3(k));
        dzL = abs(zL2(j) - zL3(k));
    
        if (dxL + dyL + dzL < min_distance)
            min_distance = dxL + dyL + dzL;
            min_xl = xL2(j);
            min_yl = yL2(j);
            min_zl = zL2(j);
        end
        
        tol = 0.0001;
        if (dxL < tol && dyL < tol && dzL < tol)
            line_intersection_points(i,1) = xL2(j);
            line_intersection_points(i,2) = yL2(j);
            line_intersection_points(i,3) = zL2(j);
            i = i + 1;
        end
    end
end
k = 2;
for i = 2:size(line_intersection_points,1)
    if sum(line_intersection_points(k,:)-line_intersection_points(k-1,:) < 0.001)
        line_intersection_points(k,:) = [];
    else
        k = k + 1;
    end
end

% visualize the intersection lines.
fig = figure;
hold on;
dl = 2;
L1 = line(xL1(1:dl:end), yL1(1:dl:end), zL1(1:dl:end), 'Color', 'b', 'LineWidth', 2);
L2 = line(xL2(1:dl:end-1000), yL2(1:dl:end-1000), zL2(1:dl:end-1000), 'Color', 'r', 'LineWidth', 2);
L3 = line(xL3(1:dl:end), yL3(1:dl:end), zL3(1:dl:end), 'Color', 'k', 'LineWidth', 2);
col = [0.3, 0.3, 0.3];
plot3(line_intersection_points(1), line_intersection_points(2), line_intersection_points(3),'o','LineWidth',2,'color',col,'MarkerSize',15);
view([-60, 20]); camlight; axis vis3d 
xlabel('[Activator] (a.u.)');
ylabel('[Inhibitor] (a.u.)');
zlabel('Active Erk (a.u.)');
FigName = 'figures/nullclines_source';
standardizePlot_3d(gcf,gca,FigName);
close(fig);

