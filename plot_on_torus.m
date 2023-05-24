%% Plot a given time series on a torus
% Note: expects a (n x 4) isothermal time series of the form (u1,u2,v1,v2)

% INPUT
% U = (u1,u2) toroidal time series (n x 2)
% V = (v1,v2) poloidal time series (n x 2)
% params = {N,q,p,c,r,a,R,cap,phi,theta,gr} (cell type)

% OUTPUT
% Plot on torus

% EXAMPLE CALL
% params = {N,q,p,c,r,a,R,cap,phi,theta,gr};
% plot_on_torus(U,V,params)

function [] = plot_on_torus(U,V,params)

% Name parameters
[~,~,~,c,r,a,R,~,~,theta,gr] = params{:};

% Font-size, for plotting
fs = 18;

% Unwrap
U = UVwrap(U, [-pi*c, pi*c]);
V = UVwrap(V, [-c*gr, c*gr]);

% Toroidal-poloidal coordinates
Phi = U./c;
Theta = [theta(V(:,1)), theta(V(:,2))];

% 3D Cartesian coordinates
X = (a + r*cos(Theta)).*cos(Phi);
Y = (R + r*cos(Theta)).*sin(Phi);
Z = r*sin(Theta);

% The parameterized torus
[Utorus, Vtorus] = meshgrid(linspace(0,2*pi), linspace(0,2*pi));
Xtorus = (a + r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R + r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);

% Plotting orbit on torus
figure (99)
surf(Xtorus, Ytorus, Ztorus,'FaceAlpha',0.3);
colormap('gray')
shading interp;
hold on
plot3(X,Y,Z)
plot3(X(1,:), Y(1,:), Z(1,:),'ok','MarkerFaceColor','k','MarkerSize',4)
plot3(X(end,:), Y(end,:), Z(end,:),'ok','MarkerFaceColor','auto','MarkerSize',4)
hold off

% Plot settings
grid on
axis equal
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)
ax = gca;
ax.FontSize = fs - 2;
set(ax,'TickLabelInterpreter','latex')
legend('','','','Start','End','Interpreter','latex','FontSize',fs-5)

end
