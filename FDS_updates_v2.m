%% Numerical NLS vortex dynamics
% On the surface of the elongated torus
clear, clc;

% Fontsize, for plotting
fs = 14;

% NLS density
mu = 3;

% Grid size (NxN)
N = 150;

% Torus parameters
a = 11;
R = 11;
r = 8;
c = sqrt(R^2 - r^2);
myalpha = r/R;

%% Solving for isothermal coordinates
% Need this to get nome
% Boundary condition
phi_0 = 0;

% Numerical integration using ode45
theta_span = [0, pi];
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);
[theta_raw, phi_raw] = ode15s(@(theta,phi) odefcn(theta,phi,a,R,r), theta_span, phi_0, options);

% Give these hosses odd symmetry 
theta_raw = [-flip(theta_raw(2:end)); theta_raw];
phi_raw = [-flip(phi_raw(2:end)); phi_raw];

% Quick plot to verify isothermal solutions look good
v0 =@(theta) 2*r*atan(sqrt((R-r)/(R+r))*tan(theta/2));

figure (1)
plot(theta_raw, real(phi_raw))
hold on
plot(theta_raw, imag(phi_raw))
plot(theta_raw, -v0(theta_raw)/c,'--k')
hold off
grid on

title('Isothermal Coordinates ODE','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$','Interpreter','latex','FontSize', fs)
ylabel('$\phi = f(\theta)$ solution','Interpreter','latex','FontSize', fs)
legend('Re$[f(\theta)]$','Im$[f(\theta)]$','$-v_0/c$','Interpreter','latex','FontSize', fs)

% Define (phi,theta) to (u,v) map
vhelper = fit(theta_raw, imag(phi_raw), 'cubicinterp');
u =@(phi) c*phi;
v =@(theta) -c*vhelper(theta);

% Defining the (u,v) to (phi,theta) map
ginv = fit(imag(phi_raw), theta_raw, 'cubicinterp');
phi =@(u) u/c;
theta =@(v) ginv(-v/c);

% Define the derivative of g_inverse
D = differentiate(ginv, imag(phi_raw));
Dginv = fit(imag(phi_raw), D, 'cubicinterp');

% Define gl and gr
gl = min(imag(phi_raw));
gr = max(imag(phi_raw));

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
w1 = (6) + 1i*(0); % positive vortex
w2 = (-6) + 1i*(0); % negative vortex

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

figure (3)
contour(Utemp,Vtemp,Z,50)
colormap hsv;
axis equal;
title('Initial Phase Contours')

% Create initial wave function
IC =@(w,w1,w2) sqrt(mu)*exp(1i*phase(w,w1,w2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w1)).^2 + (imag(w)-imag(w1)).^2))...
    .*tanh(sqrt(mu)*sqrt((real(w)-real(w2)).^2 + (imag(w)-imag(w2)).^2));
% IC =@(w,w1,w2) sqrt(mu)*exp(1i*phase(w,w1,w2))...
%     .*(tanh(0.8.*sqrt(mu)*sqrt((real(w)-real(w1)).^2 + (imag(w)-imag(w1)).^2)).^(1.5))...
%     .*(tanh(0.8.*sqrt(mu)*sqrt((real(w)-real(w2)).^2 + (imag(w)-imag(w2)).^2)).^(1.5));

% Plot density of initial condition
figure (4)
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

%% Numerical integration
% RK-4 for time
% CFL is something like dt < (dx)^2/sqrt(2) for 2D
%dt = 0.0005;
dt = 0.0025;
%dt = 0.0005;
tf = 500;
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
        disp('Plotting frame')
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
        pause(1)
    
        % Export images to folder
        if export_bool == true
            file_name = sprintf('PDE_%d.png', plot_counter);
            exportgraphics(gcf,strcat(working_dir,file_name));
        end
    end
    
    tic;
    k1 = RHS(psi, L, du, dv);
    k2 = RHS(psi + (dt/2)*k1, L, du, dv);
    k3 = RHS(psi + (dt/2)*k2, L, du, dv);
    k4 = RHS(psi + dt*k3, L, du, dv);
    psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    toc;

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

% % Wrap U,V to interval
% % IS THIS NECESSARY?
% function wrapped = UVwrap(array, interval)
%     wrapped = mod(array - interval(1), range(interval)) + interval(1);
% end

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
