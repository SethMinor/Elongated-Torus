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

E0 = -1.89;
plot_counter = 0;

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

% ode45 with events function
t0 = 0;
tf = 25000;
timespan = [t0, tf];
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11, 'Events', @EventsFcn);

% Create list of isosurface ICs (u2 = constant)
% 2 x N matrix of complex numbers (maybe N=5 or 6ish)
% Vortex 1 ICs (positive)
orbit_list(1,:) = [-3.24749, -2.49752, -2.14838, -3.427, 0.0525085, 0.0525085, -1.14749, -1.52228, -1.44749,...
	-5.99332, -4.95324, -5.04749, -6.31649, -2.94749, -2.34749, -1.74749, -2.94058,...
	2.20402, 2.15251, -0.547492, 1.02132, 1.4894, 3.65251, 7.00212, 6.67039, 5.45251, 5.15251,...
	1.85251, 2.17401, 2.40389, 1.25251, 5.26865, 4.88566, 4.56937, 5.15251,...
	8.06429, 7.01917, 6.99138, 8.10138, -8.04749, -6.27752, -6.20249, -8.23389,...
	-2.34749, -3.98657, -3.93271, -5.39858, -0.847492, -1.74749, -1.33159, -0.247492]...
    + 1i*[-0.709783, -6.42478, -7.42478, -8.92478, -3.72478, -7.81111, -4.65016, -5.52478, -6.3803,...
	2.87522, 0.175222, -2.82478, -7.32478, -0.124778, -3.72478, -6.91623, -7.62478,...
	-6.12478, -7.92478, -8.82478, -4.62478, -6.42478, -1.72557, 4.37522, 3.17522, 0.0334414, -3.81826,...
	-0.124778, -2.22478, -5.22478, -9.25687, 1.37522, -0.424778, -3.72478, -6.79009,...
	7.37522, 4.07522, -4.62478, -8.82478, 7.24806, 2.27522, -3.42478, -8.82478,...
	0.0747632, -0.124778, -3.12478, -7.92478, -3.72478, -4.85711, -6.72478, -7.77042];

% Vortex 2 ICs (negative)
orbit_list(2,:) = [0, 0,0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,...
	0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0]...
    + 1i*[7.07522, 8.57522, 6.57522, 8.07522, 7.61146, 4.07522, 7.07522, 5.87522, 5.27522,...
	7.97522, 5.27522, 0.98804, -3.12478, 8.09452, 5.384, 2.57522, 0.475222,...
	8.27522, 6.24491, 0.882667, 7.37522, 5.57522, 5.27522, 7.97522, 5.27522, 2.87522, 0.175222,...
	8.99847, 6.77522, 3.77522, 0.175222, 6.77522, 4.37522, 1.37522, -1.32478,...
	7.97522, 5.27522, -4.32478, -7.62478, 8.27522, 4.37522, -2.22478, -8.52478,...
	8.87522, 6.77522, 3.17522, -2.22478, 7.24541, 5.57522, 3.77522, 2.57522];

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

    if isempty(Ye) == 1
        fprintf('Orbit %.f, no intersections found!\n', plot_counter)
    end
    
    if isempty(Ye) == 0
        % Plot Poincare section
        % Coordinates at Poincare crossings
        Ue = Ye(:,1:N); % u-coords
        Ue = UVwrap(Ue, [-pi*c, pi*c]);
        
        Ve = Ye(:,(1+N):2*N); % v-coords
        Ve = UVwrap(Ve, [c*gl, c*gr]);

        % Check if theta1 
        
        % Conversion to toroidal-poloidal coordinates
        Phi_e = Ue./c;
        Theta_e = [theta(Ve(:,1)), theta(Ve(:,2))];
        
        % Poincare section
        figure (2)
        sgtitle("Poincare Section $(u_2 = 0)$, Energy $E_0 =$ "+E0,'Interpreter','latex')
        subplot (2,2,1)
        plot(Ue(:,1), Ue(:,2),'.','MarkerSize',3)
        xlabel('$u_1$','Interpreter','latex','FontSize',fs)
        ylabel('$u_2 = 0$','Interpreter','latex','FontSize',fs)
        xlim([-pi*c/2, pi*c/2])
        ylim([-pi*c/2, pi*c/2])
        grid on
        hold on
        
        subplot (2,2,2)
        plot(Ve(:,1), Ve(:,2),'.','MarkerSize',3)
        xlabel('$v_1$','Interpreter','latex','FontSize',fs)
        ylabel('$v_2$','Interpreter','latex','FontSize',fs)
        xlim([c*gl, c*gr])
        ylim([c*gl, c*gr])
        grid on
        hold on

        subplot (2,2,3)
        plot(Ue(:,1), Ve(:,1),'.','MarkerSize',3)
        xlabel('$u_1$','Interpreter','latex','FontSize',fs)
        ylabel('$v_1$','Interpreter','latex','FontSize',fs)
        xlim([-pi*c/2, pi*c/2])
        ylim([c*gl, c*gr])
        grid on
        hold on

        subplot (2,2,4)
        plot(Ue(:,2), Ve(:,2),'.','MarkerSize',3)
        xlabel('$u_2 = 0$','Interpreter','latex','FontSize',fs)
        ylabel('$v_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi*c/2, pi*c/2])
        ylim([c*gl, c*gr])
        grid on
        hold on

        figure (3)
        plot3(Ue(:,1), Ve(:,1), Ve(:,2),'.','MarkerSize',3)
        xlabel('$u_1$','Interpreter','latex','FontSize',fs)
        ylabel('$v_1$','Interpreter','latex','FontSize',fs)
        zlabel('$v_2$','Interpreter','latex','FontSize',fs)
        xlim([-pi*c, pi*c])
        ylim([c*gl, c*gr])
        zlim([c*gl, c*gr])
        title("Last plotted: $w_1 =$ "+u1_0+"+("+v1_0+")$i$, with "+"$w_2 =$ "+u2_0+"+("+v2_0+")$i$",...
            'Interpreter','latex','FontSize',fs)
        grid on
        hold on

        figure (5)
        plot(Ue(:,1), Ve(:,1),'.','MarkerSize',4)
        title("Poincare Section $(u_2 = 0)$, Energy $E_0=$ "+E0,'Interpreter','latex')
        xlabel('$u_1$','Interpreter','latex','FontSize',fs)
        ylabel('$v_1$','Interpreter','latex','FontSize',fs)
        xlim([-c*pi/2, c*pi/2])
        ylim([c*gl, c*gr])
        hold on

        % Export images to folder
        if export_bool == true
            file_name = sprintf('Poincare_%d.png', plot_counter);
            exportgraphics(gcf,strcat(working_dir,file_name));
        end
    end
    % Finito
end

hold off

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

  % Poincare section
  %position = V(1) + V(2); % theta1 = -theta2
  %position = U(1) - U(2); % phi1 = phi2
  %position = V(1) + V(2) + U(1) + U(2); % Hyper-plane (TRY)
  %position = V(1) - pi*r/2; %- c*gr/2; % (TRY)
  %position = V(1) + V(2) + sin(V(1) + V(2)); % Implicit line in R^3
  %position = V(2); % Looks rad but has overlap
  position = U(2);

  isterminal = 0;  % Halt integration
  % Trying direction = +/- 1 seems to help with overlapping (?)
  direction = -1;   % 0 = the zero can be approached from either direction
end
