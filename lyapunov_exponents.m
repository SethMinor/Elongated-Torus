%% Simulations on the elongated torus
clear, clc, clf;

% Font size, for plotting
fs = 10;
%'Interpreter','latex','FontSize', fs

% Parameters
% SEEMS LIKE THIS CODE MATCHES OG CODE BETTER FOR LARGER ALPHA
% BATTLE BETWEEN ISOTHERMAL AND THIS (?)
% Symplectic integrator as a possible fix (?)
a = 11;
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

%% Initial conditions for vortices
% Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
% [-pi*c,pi*c] = [-33.2475, 33.2475]
% [cgl,cgr] = [-10.5830, 10.5830]
w1_0 = (0) + 1i*(8); % positive vortex
w2_0 = (0) + 1i*(-6); % negative vortex

% Vortex charges
q1 = 1;
q2 = -1;
q = [q1 q2];   % vector of vortex charges
N = length(q); % keeping track of number of vortices

% real and imaginary parts of isothermal coords
u1_0 = real(w1_0);
v1_0 = imag(w1_0);

u2_0 = real(w2_0);
v2_0 = imag(w2_0);

%% Integrate the equations of motion
% Set total time and tolerances
t0 = 0;
tf = 500;
timespan = [t0, tf];
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

% Numerical integration using ode45
y0 = [u1_0, u2_0, v1_0, v2_0];
[t,y] = ode15s('vortex_velocity_v2',timespan, y0, options,...
    N, q, r, a, R, c, p, cap, theta, Dginv, gr);

%% Change coordinates from numerical solution
% Vortices in isothermal coordinates
U = y(:,1:N); % u-coords
V = y(:,(1+N):2*N); % v-coords

% Compute energy (before UVwrap)
[energy,classic,curve,quantum] = hamiltonian_v2(U,V,N,q,p,c,r,a,R,cap,phi,theta,gr);

% Unwrap
U = UVwrap(U, [-pi*c, pi*c]);
V = UVwrap(V, [c*gl, c*gr]);

% Full complex coordinate
W = U + 1i*V;

% Toroidal-poloidal coordinates
Phi = U./c;
Theta = [theta(V(:,1)), theta(V(:,2))];

% 3D Cartesian coordinates
X = (a+r*cos(Theta)).*cos(Phi);
Y = (R+r*cos(Theta)).*sin(Phi);
Z = r*sin(Theta);

%% Plot integrated solution
% Some nice RGB colores
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];

figure (2)

% Isothermal orbit
subplot(2,1,1)
plot(U,V,'.')
grid on
xlabel('$u = $Re$(w)$','Interpreter','latex','FontSize',fs)
ylabel('$v = $Im$(w)$','Interpreter','latex','FontSize',fs)
title('Isothermal Coordinates','Interpreter','latex','FontSize',fs)
xlim([-pi*c, pi*c])
ylim([c*gl, c*gr])

% Toriodal-poloidal orbit
subplot(2,1,2)
plotwrapped(Phi(:,1),Theta(:,1),1, [-pi, pi],[-pi, pi], 0.05, bluey)
hold on
plotwrapped(Phi(:,2),Theta(:,2),1, [-pi, pi],[-pi, pi], 0.05, orangu)
hold off
xlabel('$\phi$','Interpreter','latex','FontSize',fs)
ylabel('$\theta$','Interpreter','latex','FontSize',fs)
title('Toroidal-Poloidal Coordinates','Interpreter','latex','FontSize',fs)
xlim([-pi, pi])
ylim([-pi, pi])

% 3D Cartesian surface plot of orbits
utorus = linspace(0,2*pi);
vtorus = linspace(0,2*pi);
[Utorus,Vtorus] = meshgrid(utorus,vtorus);
Xtorus = (a+r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R+r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);

figure (3)
surf(Xtorus, Ytorus, Ztorus,'FaceAlpha',0.3);
colormap('gray')
shading interp;
hold on
plot3(X,Y,Z)
hold off
grid on
axis equal
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
zlabel('$z$','Interpreter','latex','FontSize',fs)
title('3D Cartesian (Physical Path)','Interpreter','latex','FontSize',fs)

%% Display Hamiltonian
% Compute total energy
%[energy,classic,curve,quantum] = hamiltonian_v2(U,V,N,q,p,c,r,a,R,cap,phi,theta,gr);
energy_time = linspace(t0,tf,length(energy));

figure (4);

% Plot relative difference (E(t)-E(0))/E(0)
subplot(2,1,1)
plot(energy_time,((energy-energy(1))./energy(1)))
title('Energy of solution','Interpreter','latex')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\frac{H(t)-H_0}{H_0}$','Interpreter','latex')

% Plot each energy contribution
subplot(2,1,2)
plot(energy_time,energy,'-',energy_time,classic,'--')
hold on
plot(energy_time,curve,'--',energy_time,quantum,'--')
hold off
grid on
xlabel('$t$','Interpreter','latex')
ylabel('Energy contributions','Interpreter','latex')
legend('Total, $H$','Classic','Curvature','Quantum','Interpreter','latex')

%% Compute Lyapunov exponents
% Roll through 'y_n' solution and store {r^n} values
skip_every = 1;

% Define the RHS
F =@(W) vortex_velocity_v2(0,[W(1), W(2), W(3), W(4)],0,N,q,r,a,R,c,p,cap,theta,Dginv,gr);

% Compute Jacobian at y_n
%r_list = []; % singular values
Q_n = zeros(2*N); % Initialize Q_0

% Compute {r^n}'s
tic
for n = 2:skip_every:length(y)
    % Compute J_n
    J = myjacobian(F, y(n,:));

    % Modified GS QR this guy (?)
    %[Qnew, Rnew] = QR(Jn*Qn)
    [Q_new, R] = Gram_Schmidt_QR(J);

    % Extract singular value = diags(R) (?)
    r_list(n,:) = diag(R);

    % Discard any NaNs
    if isnan(r_list(n,:)) ~= zeros(1,4)
        r_list = r_list(n-1,:);
    end

    % Next iteration
    Q_n = Q_new;
end
toc

% Find logs and remove any Inf's
log_r_list = log(r_list);
log_r_list = log_r_list(~any(isnan(log_r_list)|isinf(log_r_list),2),:);

% Compute Lyapunov spectrum
for i = 1:2*N
    Lyapunov(i) = sum(log_r_list(:,i))./tf;
end

% Plot the spectrum
figure (5)
sgtitle('Lyapunov Spectrum','Interpreter','latex','FontSize',fs+2)
subplot(4,1,1)
plot(log_r_list(:,1)./t((end-length(log_r_list)+1):end))
hold on
yline(Lyapunov(1),'--k')
hold off
grid on
xlabel('$n$','Interpreter','latex','FontSize',fs)
ylabel('$(\lambda_1)^n$','Interpreter','latex','FontSize',fs)
title("$\lambda_1=$ "+Lyapunov(1),'Interpreter','latex','FontSize',fs)

subplot(4,1,2)
plot(log_r_list(:,2)./t((end-length(log_r_list)+1):end))
hold on
yline(Lyapunov(2),'--k')
hold off
grid on
xlabel('$n$','Interpreter','latex','FontSize',fs)
ylabel('$(\lambda_2)^n$','Interpreter','latex','FontSize',fs)
title("$\lambda_2=$ "+Lyapunov(2),'Interpreter','latex','FontSize',fs)

subplot(4,1,3)
plot(log_r_list(:,3)./t((end-length(log_r_list)+1):end))
hold on
yline(Lyapunov(3),'--k')
hold off
grid on
xlabel('$n$','Interpreter','latex','FontSize',fs)
ylabel('$(\lambda_3)^n$','Interpreter','latex','FontSize',fs)
title("$\lambda_3=$ "+Lyapunov(3),'Interpreter','latex','FontSize',fs)

subplot(4,1,4)
plot(log_r_list(:,4)./t((end-length(log_r_list)+1):end))
hold on
yline(Lyapunov(4),'--k')
hold off
grid on
xlabel('$n$','Interpreter','latex','FontSize',fs)
ylabel('$(\lambda_4)^n$','Interpreter','latex','FontSize',fs)
title("$\lambda_4=$ "+Lyapunov(4),'Interpreter','latex','FontSize',fs)

%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end

% Numerical Jacobian matrix
% See (https://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf)
function J = myjacobian(func,x)
    % Set numerical derivative parameters
    N = length(x);
    F_at_x = feval(func,x);
    epsilon = 1E-12;

    % Compute numerical derivative
    xperturb = x;
    %xperturb = x + epsilon;
    J = zeros(N);
    for i = 1:N
        xperturb(i) = xperturb(i) + epsilon;
        J(:,i) = (feval(func,xperturb) - F_at_x)/epsilon;
        xperturb(i) = x(i);
    end
end

% Modified Gram-Schmidt QR, for Lyapunov exponents
% https://www.mathworks.com/matlabcentral/fileexchange/55881-gram-schmidt-orthogonalization
function [Q, R] = Gram_Schmidt_QR(X)
    % Modified Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm) 
    % which produces the same result as [Q,R]=qr(X,0)
    % Written by Mo Chen (sth4nth@gmail.com).
    [d,n] = size(X);
    m = min(d,n);
    R = zeros(m,n);
    Q = zeros(d,m);
    for i = 1:m
        v = X(:,i);
        for j = 1:i-1
            R(j,i) = Q(:,j)'*v;
            v = v-R(j,i)*Q(:,j);
        end
        R(i,i) = norm(v);
        Q(:,i) = v/R(i,i);
    end
    R(:,m+1:n) = Q'*X(:,m+1:n);
end
