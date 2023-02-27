% Contours of vortex dipole Hamiltonian
clear, clc, clf;

% Font size, for plotting
fs = 12;

% Parameters
a = 12;
R = 12;
r = 9;
c = sqrt(R^2 - r^2);
myalpha = r/R;

%% Solving for isothermal coordinates
% Boundary condition
phi_0 = 0;

% Numerical integration using ode45
theta_span = [0, pi];
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);
[theta_raw,phi_raw] = ode15s(@(theta,phi) odefcn(theta,phi,a,R,r), theta_span, phi_0, options);

% Give these hosses odd symmetry 
theta_raw = [-flip(theta_raw(2:end)); theta_raw];
phi_raw = [-flip(phi_raw(2:end)); phi_raw];

% Quick plot to verify isothermal solutions look good
v0 =@(theta) 2*r*atan(sqrt((R-r)/(R+r))*tan(theta/2));

% Numerical isothermal coordinates
figure (1)
subplot(2,1,1)
plot(theta_raw, real(phi_raw))
hold on
plot(theta_raw, imag(phi_raw))
plot(theta_raw, -v0(theta_raw)/c,'--k')
xline(pi,'-')
xline(-pi,'-')
hold off
grid on

title('Numerical Isothermal Coordinates','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel('$\phi = f(\theta)$ solution','Interpreter','latex','FontSize', fs)
legend('Re$[f(\theta)]$','Im$[f(\theta)]$','$(-v_0)/c$','Interpreter','latex','FontSize', fs)

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp'); % phi = f(theta)
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Check derivative of numerical solution
D1 = differentiate(vhelper, theta_raw);
Df = fit(theta_raw, D1, 'cubicinterp');

subplot(2,1,2)
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
plot(theta_raw, abs(imag(-1i*r./gamma(phi_raw,theta_raw)) -  Df(theta_raw)),'.-')
grid on
set(gca, 'YScale', 'log')

title("Residual, max $=$ "+max(abs(imag(-1i*r./gamma(phi_raw,theta_raw)) -  Df(theta_raw))),'Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$, poloidal coordinate','Interpreter','latex','FontSize', fs)
ylabel("RHS - $f'(\theta)$",'Interpreter','latex','FontSize', fs)
legend('Im','Interpreter','latex','FontSize', fs)

% Conformal map
w =@(phi,theta) u(phi) + 1i*v(theta);

% Define gl and gr
gl = min(imag(phi_raw));
gr = max(imag(phi_raw));

% Defining the (u,v) to (phi,theta) map
ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
phi =@(u) u/c;
theta =@(v) ginv(-UVwrap(v, [-c*gr, c*gr])/c);

% Define the derivative of g_inverse
D2 = differentiate(ginv, imag(phi_raw));
Dginv = fit(imag(phi_raw), D2, 'cubicinterp');

%% Other parameters
% Gamma quantity
gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);

% Local scale factor
lambda =@(phi,theta) gamma(phi,theta)./c;

% Jacobi theta function parameters
cap = 20;     % truncation error
p = exp(-gr); % nome (for periodicity)

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

%% Find contours of constant energy
% % Hamiltonian as anonymous function
% E =@(phi1, phi2, theta1, theta2) hamiltonian_contour(u(phi1), u(phi2), v(theta1), v(theta2),...
%     N,q,p,c,r,a,R,cap,phi,theta,gr);
% 
% % Fix Hamiltonian in two variables
% phi1_const = 0;
% phi2_const = 1;
% E_contour =@(Theta1, Theta2) E(phi1_const, phi2_const, Theta1, Theta2);
% 
% % Plot contours of Hamiltonian
% mesh_density = 150; % Density for contour plot
% theta1 = linspace(-pi, pi, mesh_density);
% theta2 = linspace(-pi, pi, mesh_density);
% [Theta1, Theta2] = meshgrid(theta1, theta2);
% Z = E_contour(Theta1, Theta2);
% Z = reshape(Z, [mesh_density, mesh_density]); % Reshape into matrix
% 
% % Remove Infs from Z
% % Just set to NaN
% Z(isinf(Z)) = NaN;
% 
% % Normalize Z
% Z = Z./max(abs(Z),[],'all');
% 
% figure (2)
% sgtitle("Parameters: $a=$ "+a+", $R=$ "+R+", $r=$ "+r, 'interpreter', 'latex')
% 
% subplot(2,1,1)
% surf(Theta1, Theta2, Z)
% colorbar
% colormap parula;
% clim([-0.12 0])
% shading flat
% grid on
% title("Hamiltonian for $\phi_1=$ "+phi1_const+", $\phi_2=$ "+phi2_const, 'interpreter', 'latex')
% xlabel('$\theta_1$', 'interpreter', 'latex')
% ylabel('$\theta_1$', 'interpreter', 'latex')
% zlabel('Normalized Energy, $E/E_{max}$', 'interpreter', 'latex')
% 
% subplot(2,1,2)
% contour(Theta1, Theta2, Z, 2*mesh_density)
% colormap parula;
% clim([-0.12 0])
% colorbar
% grid on
% title("Contours of Hamiltonian for $\phi_1=$ "+phi1_const+", $\phi_2=$ "+phi2_const, 'interpreter', 'latex')
% xlabel('$\theta_1$', 'interpreter', 'latex')
% ylabel('$\theta_1$', 'interpreter', 'latex')
% zlabel('Normalized Energy, $E/E_{max}$', 'interpreter', 'latex')

%% Isosurfaces for circular torus
% Energy levels to check
%E0 = [-300, -200, - 10, -60, 0, 60, 175, 330, 510, 800];
E0 = linspace(-5, 5, 10);

% Circular torus gamma in isothermal
gamma_c =@(v) R + r*cos(2*atan(sqrt((R+r)/(R-r))*tan(v/(2*r))));

% Define fcn for implicit surface
f_E =@(v1,v2,d) log(abs( jacobitheta1((d + 1i*(v1-v2))./(2*c),p,cap) )) + ...
    log(abs( jacobitheta1(-(d + 1i*(v1-v2))./(2*c),p,cap) )) + ...
    log(gamma_c(v1).*gamma_c(v2)./(c^2)) + ...
    (1/(2*c^2*gr))*(d.^2);

% Solving this BL
mesh_density = 0.3;

d = -2*pi*c:mesh_density:(2*pi*c);
v1 = (-pi*r):mesh_density:(pi*r);
v2 = v1;

[V1,V2,D] = meshgrid(v1, v2, d);
V = f_E(V1,V2,D);

% Plot it!
figure (2)
for i = 1:length(E0)
    isosurface(V1, V2, D, V, E0(i));
    alpha(0.8)
    hold on
end

grid on
%xlim([0, 2*pi*c])
%ylim([-pi*r, pi*r])
%zlim([-pi*r, pi*r])

title("$H(u_1, u_2, v_1, v_2) = E_0$,  for  $(a,R,r)=($"+a+","+R+","+r+")", ...
    'interpreter', 'latex','FontSize',fs+2)
xlabel('$v_1$', 'interpreter', 'latex','FontSize',fs+2)
ylabel('$v_2$', 'interpreter', 'latex','FontSize',fs)
zlabel('$\delta = u_1 - u_2$', 'interpreter', 'latex','FontSize',fs+2)
legend("$E_0=$ "+E0(1), 'interpreter', 'latex')
% legend("$E_0=$ "+E0(1), "$E_0=$ "+E0(2),"$E_0=$ "+E0(3), "$E_0=$ "+E0(4),...
%     "$E_0=$ "+E0(5), "$E_0=$ "+E0(6), "$E_0=$ "+E0(7), "$E_0=$ "+E0(8),...
%     "$E_0=$ "+E0(9), "$E_0=$ "+E0(10), 'interpreter', 'latex')
%camlight
%lighting gouraud;

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Hamiltonian for contours
function [energy,classic,curve,quantum] = hamiltonian_contour(u1,u2,v1,v2,N,q,p,c,r,a,R,cap,phi,theta,gr)
    U = [u1, u2];
    V = [v1, v2];
    W = U + 1i*V;
    
    % Gamma
    gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
    
    % Local scale factor
    L =@(u,v) gamma(phi(u),theta(v))/c;
    
    % Initialize energy contributions
    classic = 0; curve = 0; quantum = 0;
    
    % Loop over each vortex pair
    for n = 1:N
        % Curvature contribution
        curve = curve + log(L(U(:,n),V(:,n)));
    
        for m = 1:N
            % quantum contribution
            quantum = quantum + (1/(2*c^2*gr)).*q(n).*q(m).*U(:,n).*U(:,m);
    
            if (m ~= n)
                % classical contribution
                thetus = (W(:,n) - W(:,m))./(2*c);
                classic = classic -...
                    (q(n).*q(m).*log(abs(jacobitheta1(thetus,p,cap))));
            end
        end
    end
    energy = classic + curve + quantum;
end
