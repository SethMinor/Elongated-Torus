%% Numerical NLS vortex dynamics
% On the surface of the elongated torus
clear, clc;

% Fontsize, for plotting
fs = 14;

% NLS density
mu = 4;

% Grid size (NxN)
N = 250;

% Torus parameters
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
%theta_span = [0, pi+-0.01];
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
plot(theta_raw, imag(-1i*r./gamma(phi_raw,theta_raw)) -  Df(theta_raw),'.-')
grid on

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

%% Initial conditions for wave function
% Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
w1 = (0) + 1i*(3); % positive vortex
w2 = (0) + 1i*(-3); % negative vortex

% Complex flow potential
F =@(w,w1,w2) log(jacobitheta1((w-w1)./(2*c),p,cap)./jacobitheta1((w-w2)./(2*c),p,cap))...
    - (real(w1-w2)./(2*c^2*gr))*w;

% Phase of the flow field, Phi=Im(F)
phase =@(w,w1,w2) imag(F(w,w1,w2));

% Create a contour plot of the phase
du = 2*pi*c/N;
dv = 2*c*gr/N;
Ugrid = linspace(-c*pi,c*pi - du,N);
Vgrid = linspace(c*gl,c*gr - dv,N);
[Utemp,Vtemp] = meshgrid(Ugrid,Vgrid);
Z = phase(Utemp+1i*Vtemp,w1,w2);

figure (2)
contour(Utemp,Vtemp,Z,50)
colormap hsv;
axis equal;
title('Initial Phase Contours')

% Create initial wave function (PDE)
IC =@(w,w1,w2) sqrt(mu)*exp(1i*phase(w,w1,w2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w1)).^2 + (imag(w)-imag(w1)).^2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w2)).^2 + (imag(w)-imag(w2)).^2));

% Plot density of initial condition
figure (3)
psi_0 = IC(Utemp+1i*Vtemp,w1,w2);
Z = conj(psi_0).*psi_0;
surf(Utemp,Vtemp,Z)
shading interp;
colormap gray;
axis equal;
camlight
title('Initial Density')

% For plotting on el toroos:
% https://www.mathworks.com/help/matlab/visualize/representing-a-matrix-as-a-surface.html

%% Numerical integration (ODE)
% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N_vortices = length(q); % keeping track of number of vortices

% Real and imaginary parts of ICs
u1 = real(w1);
v1 = imag(w1);

u2 = real(w2);
v2 = imag(w2);

% Set equations of motion
tf = 800;
timespan = [0, tf];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Numerical integration using ode45
y0 = [u1, u2, v1, v2];
[t_ode, y_ode] = ode15s('vortex_velocity_v2',timespan, y0, options,...
    N_vortices, q, r, a, R, c, p, cap, theta, Dginv, gr);

% Vortices in isothermal coordinates
U = y_ode(:,1:N_vortices); % u-coords
V = y_ode(:,(1+N_vortices):2*N_vortices); % v-coords

% Compute energy before UVwrap
[energy,classic,curve,quantum] = hamiltonian_v2(U,V,N_vortices,...
    q,p,c,r,a,R,cap,phi,theta,gr);

% UVwrap time
U = UVwrap(U, [-pi*c, pi*c]);
V = UVwrap(V, [c*gl, c*gr]);

% Some nice RGB colores
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];
figure (4)

% Physical trajectory in isothermal coordinates
subplot(3,1,1)
plot(U,V,'.')
grid on
xlabel('$u = $Re$(w)$','Interpreter','latex','FontSize',fs)
ylabel('$v = $Im$(w)$','Interpreter','latex','FontSize',fs)
title('Isothermal Coordinates','Interpreter','latex','FontSize',fs)
xlim([-pi*c, pi*c])
ylim([c*gl, c*gr])

% Energy conservation
% [energy,classic,curve,quantum] = hamiltonian_v2(U,V,N_vortices,...
%     q,p,c,r,a,R,cap,phi,theta,gr);
energy_time = linspace(0,tf,length(energy));

% Plot relative difference (E(t)-E(0))/E(0)
subplot(3,1,2)
plot(energy_time,((energy-energy(1))./energy(1)))
title('Energy of solution','Interpreter','latex')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\frac{H(t)-H_0}{H_0}$','Interpreter','latex')

% Plot each energy contribution
subplot(3,1,3)
plot(energy_time,energy,'-',energy_time,classic,'--')
hold on
plot(energy_time,curve,'--',energy_time,quantum,'--')
hold off
grid on
xlabel('$t$','Interpreter','latex')
ylabel('Energy contributions','Interpreter','latex')
legend('Total, $H$','Classic','Curvature','Quantum','Interpreter','latex')

%% Numerical integration (PDE)
% RK-4 for time
% CFL is something like dt < (dx)^2/sqrt(2) for 2D
%dt = 0.0005;
dt = 0.001;
% dt = 0.005 % Largest dt found to still work
N_time = floor(tf/dt);

% Local scale factor
figure (5)
subplot(1,2,1)
[Phi_temp, Theta_temp] = meshgrid(phi(Ugrid),theta(Vgrid)');
surf(Phi_temp, Theta_temp, lambda(Phi_temp,Theta_temp))
title('\Lambda(\phi,\theta)')
xlabel('\phi')
ylabel('\theta')
view([0,90])
shading interp;
colormap default;

subplot(1,2,2)
surf(Utemp, Vtemp, lambda(Phi_temp,Theta_temp))
title('\Lambda(u,v)')
xlabel('u')
ylabel('v')
view([0,90])
shading interp;
colormap default;

% Reshape initial condition array
seed = psi_0;

% Define spatially-varying Lambda as matrix
% Multiply by pre-computed 1/L instead of dividing for performance (?)
L = lambda(Phi_temp,Theta_temp);

% Want to export images?
export_bool = true;
working_dir = 'C:\Users\sminor2848\Downloads\Elongated-Torus-main\Elongated-Torus-main\pics\';

% RK-4 for-loop
t = 0; % Initialize time
psi = seed; % Initialize wave function

% Confirm IC
disp('Continue? Press ye olde enter key...')
pause;
plot_counter = 0;

% Initialize mass history
mass = [];

for i = 0:N_time
    % Plot every 0.5ish seconds
    if mod(i,1000) == 0
        plot_counter = plot_counter + 1;
        figure (6)
        %disp('Plotting frame')
        % Plot time step
        density = conj(psi).*psi;
        % Possibly a different mass formula due to curvature?
        mass(plot_counter) = sum(sum(conj(psi).*psi));
        surf(Utemp,Vtemp,density)
        shading interp;
        colormap gray;
        axis equal;
        view(0,90)
        colorbar
        %camlight
        title("$t=$ "+t,'Interpreter','latex','FontSize',fs)
        xlim([-pi*c,pi*c])
        ylim([c*gl,c*gr])
        %pause(1)
        
        % Add ODE orbit
        hold on
        % Get index of point at the nearest time
        [~,i_nearest] = (min(abs(t_ode - t)));
        plot3(U(1:i_nearest,:), V(1:i_nearest,:), 0*U(1:i_nearest,:)+mu+1,...
            '.','MarkerSize',4)
        hold off
    
        % Export images to folder
        if export_bool == true
            file_name = sprintf('PDE_%d.png', plot_counter);
            exportgraphics(gcf,strcat(working_dir,file_name));
        end
    end
    
    %tic;
    k1 = RHS(psi, L, du, dv);
    k2 = RHS(psi + (dt/2)*k1, L, du, dv);
    k3 = RHS(psi + (dt/2)*k2, L, du, dv);
    k4 = RHS(psi + dt*k3, L, du, dv);
    psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    %toc;

    % Update the time
    t = t + dt;
end

% Display mass conservation
figure (7)
plot((mass-mass(1))/mass(1), '--xr')
title('Mass Conservation')
xlabel("n, where t = (n-1)"+1000*dt)
ylabel('Relative Error')

%% Helper functions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Return the RHS
function F_of_psi = RHS(psi, L, du, dv)
    %F_of_psi = 1i.*(2*del2(psi, du, dv)./L - (abs(psi).^2).*psi);

    % Periodic BCs w/ circshift?
    % Advice: don't use circshift for performance (?) home-made instead (?)
    % Check 639 code
    % Define 1/du^2 or dv^2 as variable to speed up performance (?)
    % (pre-computing)
    % Matlab 'profiler' performance GUI - google, yo!
    % Rewrite expression to have only one multiplication if possible
    F_of_psi = (1i/2)*(-2*(1/(du^2) + 1/(dv^2))*psi ...
    + (1/(du^2))*( circshift(psi,1,2) + circshift(psi,-1,2) )...
    + (1/(dv^2))*( circshift(psi,1,1) + circshift(psi,-1,1) ))./(L.^2) ...
    - 1i*(conj(psi).*psi).*psi;
    % Maybe try psi.^2 (?)
end
