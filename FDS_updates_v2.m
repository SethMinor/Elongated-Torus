%% Numerical NLS vortex dynamics
% On the surface of the elongated torus
clear, clc;

% Fontsize, for plotting
fs = 14;

% NLS density
mu = 3;

% Grid size (NxN)
N = 50;

% Torus parameters
a = 11;
R = 11;
r = 3;
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

% Plot to verify that 'v' function looks good
figure (2)
plot(theta_raw, v(theta_raw))
hold on
plot(theta_raw, v0(theta_raw),'--k')
hold off
grid on

title('Computed $v$ versus known $v_0$','Interpreter','latex','FontSize', fs+2)
xlabel('$\theta$','Interpreter','latex','FontSize', fs)
ylabel('$v(\theta)$','Interpreter','latex','FontSize', fs)
legend('$v$','$v_0$','Interpreter','latex','FontSize', fs)

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
w1 = (0) + 1i*(4); % positive vortex
w2 = (0) + 1i*(-4); % negative vortex

% Complex flow potential
F =@(w,w1,w2) log(jacobitheta1((w-w1)./(2*c),p,cap)./jacobitheta1((w-w2)./(2*c),p,cap))...
    - (real(w1-w2)./(2*c^2*gr))*w;

% Phase of the flow field, Phi=Im(F)
phase =@(w,w1,w2) imag(F(w,w1,w2));

% Create a contour plot of the phase
Ugrid = linspace(-c*pi,c*pi,N);
Vgrid = linspace(c*gl,c*gr,N);
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
dt = 0.01;
tf = 100;
N_time = floor(tf/dt);

% Local scale factor
figure (5)
subplot(1,2,1)
[Phi_temp, Theta_temp] = meshgrid(phi(Ugrid),theta(Vgrid)');
surf(Phi_temp, Theta_temp, lambda(Phi_temp,Theta_temp))
title('Local scale factor, \Lambda')
xlabel('\phi')
ylabel('\theta')
view([0,90])
shading interp;
colormap default;

subplot(1,2,2)
surf(Utemp, Vtemp, lambda(Phi_temp,Theta_temp))
title('\Lambda (isothermal)')
xlabel('u')
ylabel('v')
view([0,90])
shading interp;
colormap default;

% Reshape initial condition array
seed = psi_0;

% Wrap spatially-varying Lambda into vector
L = lambda(Phi_temp,Theta_temp);

% Confirm IC
disp('Continue? Press ye olde enter key...')
pause;

% Want to export images?
export_bool = true;
working_dir = 'C:\Users\sminor2848\Downloads\Elongated-Torus-main\Elongated-Torus-main\pics\';

% RK-4 for-loop
figure (6)
t = 0; % Initialize time
psi = seed; % Initialize wave function

% Discretization parameters
hu = 2*pi*c/N;
hv = 2*c*gr/N;

for i = 0:N_time
    % Plot time step
    density = conj(psi).*psi;
    mass = sum(density,'all');
    surf(Utemp,Vtemp,density)
    shading interp;
    colormap gray;
    axis equal;
    view(0,90)
    colorbar
    %camlight
    title("$t=$ "+t+", $M=$ "+mass,'Interpreter','latex','FontSize',fs)
    xlim([-pi*c,pi*c])
    ylim([c*gl,c*gr])
    %pause(0.01)
    %contour(Utemp,Vtemp,phase(psi,w1,w2),10)
    %colormap hsv;
    %axis equal;

    % Export images to folder
    if export_bool == true
        file_name = sprintf('PDE_%d.png', i);
        exportgraphics(gcf,strcat(working_dir,file_name))
    end
    
    tic;
    k1 = RHS(psi, L, hu, hv);
    k2 = RHS(psi + (dt/2)*k1, L, hu, hv);
    k3 = RHS(psi + (dt/2)*k2, L, hu, hv);
    k4 = RHS(psi + dt*k3, L, hu, hv);
    psi = psi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    toc;
    %psi = psi + dt*RHS(psi, L, hu, hv);

    % Enforce periodic BCs?
    psi(N,:) = psi(1,:);
    psi(:,N) = psi(:,1);

    % Update the time
    t = t + dt;
end


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
function F_of_psi = RHS(psi, L, hu, hv)
    F_of_psi = 1i.*(4*del2(psi, hu, hv)./L - (abs(psi).^2).*psi);
end
