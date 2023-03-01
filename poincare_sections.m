% Poincare sections for elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;

% Parameters
a = 13;
R = 11;
r = 3;
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
%theta =@(v) ginv(-v/c);
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

%% Giant for-loop over ICs
% Want to export images?
export_bool = false;
working_dir = 'C:\Users\sminor2848\Downloads\Elongated-Torus-main\Elongated-Torus-main\pics\';

plot_counter = 0;

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

% ode45 with events function
t0 = 0;
tf = 10000;
timespan = [t0, tf];
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11, 'Events', @EventsFcn);

% Create list of isosurface ICs (u2 = constant)
% 2 x N matrix of complex numbers (maybe N=5 or 6ish)
%u2_const = 0;
% Vortex 1 ICs (positive)
orbit_list(1,:) = [3.5, -3.5 0, 0, 2.24749, 2.25251, 3.01938, -3.01938, -10.0475, -8.4, -11.2, -12.9, -13.1708, 9.9, 13.1, 7.25251, 12.8809, 9.65251, 9.95251, 9.65251, 11.8594, -11.8594]...
    + 1i*[0, 0, 2, -2, -1.49626, 1.64215, 0.924778, 0.924778, -1.82478, 1.6, 1.6, 0.8, 0.175222, -1.7, -0.1, -1.02478, -0.724778, 3.77522, -8.22478, 9, -0.024778, -0.024778];
 % Vortex 2 ICs (negative)
orbit_list(2,:) = [0, 0, 0, 0, 0, 0, 0, 0, -10, -10, -10, -10, -10, 10, 10, 10, 10, 10, 10, 10, 15, -15]...
    + 1i*[0, 0, -2, 2, 1.57522, -1.42478, -0.924778, -0.924778, 1.77439, -1.6, -1.6, -0.8, -0.124778, 1.9, 0.1, 0.962687, 0.775222, 8.61893, -3.46464, 4.075, -0.024778, -0.024778];

% REMEMBER that shown surfaces are just slices of full 3D volume
% So ICs on the shown surface may 4D rotate into u2 =/= 0

% Periodic orbits for u1=u2, v1=-v2

for IC_number = 1:length(orbit_list)
    plot_counter = plot_counter + 1;
    % Initial conditions for vortices
    % Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
    w1_0 = orbit_list(1, IC_number); % positive vortex
    w2_0 = orbit_list(2, IC_number); % negative vortex
    
    % real and imaginary parts of isothermal coords
    u1_0 = real(w1_0);
    v1_0 = imag(w1_0);
    
    u2_0 = real(w2_0);
    v2_0 = imag(w2_0);

    % Compute Poincare sections
    % Numerical integration using ode45 or ode15s
    y0 = [u1_0, u2_0, v1_0, v2_0];
    [T,Y,Te,Ye,Ie] = ode15s(@(t,y) vortex_velocity_v2(0,y,0,N,q,r,a,R,c,p,cap,theta,Dginv,gr),...
        timespan, y0, options);
    
    if isempty(Ye) == 0
        % Plot Poincare section
        % Coordinates at Poincare crossings
        Ue = Ye(:,1:N); % u-coords
        Ue = UVwrap(Ue, [-pi*c, pi*c]);
        
        Ve = Ye(:,(1+N):2*N); % v-coords
        Ve = UVwrap(Ve, [c*gl, c*gr]);
        
        % Conversion to toroidal-poloidal coordinates
        Phi_e = Ue./c;
        Theta_e = [theta(Ve(:,1)), theta(Ve(:,2))];
        
        % Poincare section
        figure (2)
        sgtitle("Poincare Section $(\theta_1 = -\theta_2)$, Energy $E_0 = -3$",'Interpreter','latex')
        subplot (2,2,1)
        plot(Phi_e(:,1), Phi_e(:,2),'.','MarkerSize',3)
        xlabel('$\phi_1$','Interpreter','latex','FontSize',fs)
        ylabel('$\phi_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        grid on
        hold on
        
        subplot (2,2,2)
        plot(Theta_e(:,1), Theta_e(:,2),'.','MarkerSize',3)
        xlabel('$\theta_1$','Interpreter','latex','FontSize',fs)
        ylabel('$\theta_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        grid on
        hold on

        subplot (2,2,3)
        plot(Phi_e(:,1), Theta_e(:,1),'.','MarkerSize',3)
        xlabel('$\phi_1$','Interpreter','latex','FontSize',fs)
        ylabel('$\theta_1 = -\theta_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        hold on

        subplot (2,2,4)
        plot(Phi_e(:,2), Theta_e(:,2),'.','MarkerSize',3)
        xlabel('$\phi_2$','Interpreter','latex','FontSize',fs)
        ylabel('$\theta_2 = -\theta_1$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        hold on

        figure (3)
        plot3(Phi_e(:,1), Phi_e(:,2), Theta_e(:,1),'.','MarkerSize',3)
        xlabel('$\phi_1$','Interpreter','latex','FontSize',fs)
        ylabel('$\phi_2$','Interpreter','latex','FontSize',fs)
        zlabel('$\theta_1 = -\theta_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        zlim([-pi, pi])
        title("Last plotted: $w_1 =$ "+u1_0+"+("+v1_0+")$i$, with "+"$w_2 =$ "+u2_0+"+("+v2_0+")$i$",...
            'Interpreter','latex','FontSize',fs)
        grid on
        hold on

        figure (5)
        plot(Phi_e(:,1), Theta_e(:,1),'.k','MarkerSize',2)
        title("Poincare Section $(\theta_1 = -\theta_2)$, Energy $E_0 = -3$",'Interpreter','latex')
        xlabel('$\phi_1$','Interpreter','latex','FontSize',fs)
        ylabel('$\theta_1 = -\theta_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi, pi])
        ylim([-pi, pi])
        hold on

        % Export images to folder
        if export_bool == true
            file_name = sprintf('Poincare_%d.png', plot_counter);
            exportgraphics(gcf,strcat(working_dir,file_name));
        end
    end
    % Finito
end

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Poincare section slice
% UPDATE THIS BL WHEN R,r CHANGE (!!!)
function [position,isterminal,direction] = EventsFcn(~,y)
  % Torus parameters
  R = 11;
  r = 3;
  c = sqrt(R^2-r^2);
  gr = 1.0345;
  % Isothermal coordinates
  U = y(1:2);
  U = UVwrap(U, [-pi*c, pi*c]);
  V = y(3:4);
  V = UVwrap(V, [-c*gr, c*gr]);
  % Toroidal-poloidal cooridnates
  %Phi = U./c;
  % Poincare section
  position = V(1) + V(2); % theta1 = -theta2
  %position = U(1) - U(2); % phi1 = phi2
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end
