% Poincare sections for elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 18;

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

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp'); % phi = f(theta)
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Check derivative of numerical solution
D1 = differentiate(vhelper, theta_raw);
Df = fit(theta_raw, D1, 'cubicinterp');

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
tf = 30000;
timespan = [t0, tf];
options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11, 'Events', @EventsFcn);

% Create list of isosurface ICs (u2 = constant)
% 2 x N matrix of complex numbers (maybe N=5 or 6ish)
% Vortex 1 ICs (positive)
orbit_list(1,:) = [0.352508, 0.352508, 2.45251, 2.45251, 2.25251, -2.25251, 2.15251, -2.15251,...
	1.25251, 0.352508, -5.32283, -2.67955, -5.23302, 5.32368, 2.45251, 5.19891,...
	3.35251, 2.91812, 2.45251, 0.0525085, 0.0525085, -0.547492, -2.64749, -3.12493,...
	2.45251, 3.82508, -0.247492, 0.0525085, -2.34749, -3.84749, 3.05251, 0.352508, -0.847492,...
	1.25251, -1.44749, -3.84749, 5.24173, 4.89834, 5.01941, 5.39979, -5.14663, -4.89729, -4.94494, -5.28658,...
    3.8475]...
    + 1i*[8.59535, 2.57522, -4.17478, -4.17478, -4.97819, -4.97819, -5.72478, -5.72478,...
	2.87522, 2.41657, 7.37522, -2.52478, -8.22478, 7.37522, -2.52478, -7.62478,...
	4.5782, -0.124778, -3.32573, 5.72388, 1.26486, -1.92478, 7.07522, 0.475222,...
	7.07522, -4.02478, 5.87522, -2.89081, 6.47522, -4.32478, 7.37522, 6.94037, 5.87522,...
	-2.62336, -2.68035, -3.72478, 8.27522, 8.27522, -6.72478, -7.92478, 8.27522, 7.67522, -5.82478, -7.62478,...
    -3.7248];

% Vortex 2 ICs (negative)
orbit_list(2,:) = [0,0,0, 0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0, 0,0,0, 0,0,0, ...
	0,0,0, 0,0,0, 0,0,0,0, 0,0,0,0,...
    0, 0,0,0, 0,0]...
    + 1i*[3.17522, -1.4753, -9.12478, -9.12478, -1.22478, -1.22478, 1.67385, 1.67385,...
	-0.955982, -1.62478, 8.27522, 0.175222, -7.02478, 8.27522, 0.425142, -6.72478,...
	1.97522, -2.52478, -7.92478, 1.07522, -2.82478, -6.84328, 2.8801, -1.32478,...
	2.75726, -7.32478, 1.19285, -8.22478, 2.29367, -7.77856, 3.36339, 1.97522,...
	1.27455, -7.62478, -7.62478, -6.75811, 7.07522, 6.17522, -8.82478, -8.82478, 6.77522, 5.87522, -7.32478, -9.12478,...
    -6.7581];

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
        figure (1)
%         sgtitle("Poincare Section $(u_2 = 0)$, Energy $E_0 =$ "+E0,'Interpreter','latex')
%         subplot (2,2,1)
%         plot(Ue(:,1), Ue(:,2),'.','MarkerSize',3)
%         xlabel('$u_1$','Interpreter','latex','FontSize',fs)
%         ylabel('$u_2 = 0$','Interpreter','latex','FontSize',fs)
%         xlim([-pi*c/2, pi*c/2])
%         ylim([-pi*c/2, pi*c/2])
%         grid on
%         hold on
        
%         subplot (2,2,2)
%         plot(Ve(:,1), Ve(:,2),'.','MarkerSize',3)
%         xlabel('$v_1$','Interpreter','latex','FontSize',fs)
%         ylabel('$v_2$','Interpreter','latex','FontSize',fs)
%         xlim([c*gl, c*gr])
%         ylim([c*gl, c*gr])
%         grid on
%         hold on

        %subplot (2,2,3)
        plot(Ue(:,1), Ve(:,1),'.','MarkerSize',5)
        xlabel('$u_1$','Interpreter','latex','FontSize',fs)
        ylabel('$v_1$','Interpreter','latex','FontSize',fs)
        xlim([-pi*c/4, pi*c/4])
        ylim([c*gl, c*gr])
        ax = gca;
        ax.FontSize = fs - 3;
        %grid on
        hold on

%         subplot (2,2,4)
%         plot(Ue(:,2), Ve(:,2),'.','MarkerSize',3)
%         xlabel('$u_2 = 0$','Interpreter','latex','FontSize',fs)
%         ylabel('$v_2$','Interpreter','latex','FontSize',fs)
%         xlim([-pi*c/2, pi*c/2])
%         ylim([c*gl, c*gr])
%         grid on
%         hold on

%         figure (3)
%         plot3(Ue(:,1), Ve(:,1), Ve(:,2),'.','MarkerSize',3)
%         xlabel('$u_1$','Interpreter','latex','FontSize',fs)
%         ylabel('$v_1$','Interpreter','latex','FontSize',fs)
%         zlabel('$v_2$','Interpreter','latex','FontSize',fs)
%         xlim([-pi*c, pi*c])
%         ylim([c*gl, c*gr])
%         zlim([c*gl, c*gr])
%         title("Last plotted: $w_1 =$ "+u1_0+"+("+v1_0+")$i$, with "+"$w_2 =$ "+u2_0+"+("+v2_0+")$i$",...
%             'Interpreter','latex','FontSize',fs)
%         grid on
%         hold on
% 
%         figure (5)
%         plot(Ue(:,1), Ve(:,1),'.','MarkerSize',4)
%         title("Poincare Section $(u_2 = 0)$, Energy $E_0=$ "+E0,'Interpreter','latex')
%         xlabel('$u_1$','Interpreter','latex','FontSize',fs)
%         ylabel('$v_1$','Interpreter','latex','FontSize',fs)
%         xlim([-c*pi/2, c*pi/2])
%         ylim([c*gl, c*gr])
%         hold on
        
        fprintf('Just printed orbit %.f\n',IC_number)

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
  position = U(2);

  isterminal = 0;  % Halt integration
  % Trying direction = +/- 1 seems to help with overlapping (?)
  direction = 1;   % 0 = the zero can be approached from either direction
end
