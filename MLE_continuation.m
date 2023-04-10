%% Branching for MLE plot
clear, clc;%, clf;

% Font size, for plotting
fs = 18;

% Compute Lyapunov exponents as well?
Lyap_bool = false;

% Initial parameters
a = 13.95;
a0 = a;
R = 11;
r = 2.64;
c = sqrt(R^2 - r^2);
myalpha = r/R;

% Energy level to check
E0 = -9.87844;

% Initial seeds for [v1,v2]
v_seed = [-10, -9.8];

% Set u2 = constant = u1
u2constant = -30;

%% Select a branch of ICs for testing
% Continuation for-loop
da = 0.25;
counter = 1;

for delta_a = 0:da:2.5 % (a will go from a0 to a0 - last number)
    % New semi-major axis a
    a = a0 - delta_a;

    % Initial seed
    if counter == 1
        V = v_seed;
    else
        V = v_star;
    end
    
    % --- Solving for isothermal coordinates ---
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
    
    % --- Other parameters ---
    % Gamma quantity
    gamma =@(phi,theta) sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
    
    % Local scale factor
    lambda =@(phi,theta) gamma(phi,theta)./c;
    
    % Jacobi theta function parameters
    cap = 20;     % truncation error
    p = exp(-gr); % nome (for periodicity)

    % Define implicit function f = H(w1,w2)-E0
    f =@(u1,v1,v2,E0) log(abs( jacobitheta1(((u1 - u2constant) + 1i*(v1-v2))./(2*c),p,cap) )) + ...
        log(abs( jacobitheta1(-((u1 - u2constant) + 1i*(v1-v2))./(2*c),p,cap) )) + ...
        log(gamma(u1./c, theta(v1)).*gamma(u2constant./c, theta(v2))./(c^2)) + ...
        (1/(2*c^2*gr))*((u1 - u2constant).^2) - E0;
    
    % Giving f one parameter Z = [u1,v1,v2,E0]
    f_oneparam =@(V) f(u2constant, V(1), V(2), E0);
    
    % Optimizer to solve for (v1,v2) values given seed
    f_options = optimset('Display','off','Algorithm','levenberg-marquardt');
    v_star = fsolve(f_oneparam, v_seed, f_options);

    % Update list of (a,v1,v2)
    v_list(counter,1:3) = [a,v_star];

    % Compute Lyapunov exponents of branch using continuous QR?
    if Lyap_bool == true
        % --- Initial conditions for vortices ---
        % Initial vortex positions in [-pi*c,pi*c]x[cgl,cgr]
        w1_0 = (u2constant) + 1i*(v_star(1)); % positive vortex
        w2_0 = (u2constant) + 1i*(v_star(2)); % negative vortex
        
        % Vortex charges
        q = [1, -1];   % vector of vortex charges
        N = length(q); % keeping track of number of vortices
        
        % Real and imaginary parts of isothermal coords
        u1_0 = real(w1_0);
        v1_0 = imag(w1_0);
        
        u2_0 = real(w2_0);
        v2_0 = imag(w2_0);

        % Integrate the equations of motion
        % Set total time and tolerances
        t0 = 0;
        tf = 500;
        timespan = [t0, tf];
        options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
        
        % Define the RHS
        F =@(W) vortex_velocity_v2(0,[W(1), W(2), W(3), W(4)],0,N,q,r,a,...
            R,c,p,cap,theta,Dginv,gr);
        
        % Numerical integration using ode45
        y0 = [u1_0, u2_0, v1_0, v2_0]'; % ODE ICs
        Q0 = eye(2*N); % Lyapunov, Q ICs
        p0 = ones(2*N,1); % Lyapunov, rho_i ICs
        Y0 = [y0; reshape(Q0,[(2*N)^2,1]); p0];
        
        tic
        [t,y] = ode45('vortex_velocity_lyapunov',timespan, Y0, options,...
            N, q, r, a, R, c, p, cap, theta, Dginv, gr, F);
        toc

        % Roll through 'y_n' solution and store {rho_i} values
        % [rho_1, ..., rho_4] list = y(:, 1 + 2N + (2N)^2 : end)
        rho_list = y(:, 1+2*N+(2*N)^2 : end);
        
        % Compute the full spectrum
        N_lyapunov = 100; % number to back-average
        Lyapunov = zeros(1,2*N);
        for i = 1:2*N
            Lyapunov(i) = mean(rho_list(end-N_lyapunov:end,i)./...
                t(end-N_lyapunov:end));
        end

        % Update Lyapunov list of (a, MLE)
        MLE_list(counter,1:2) = [a, Lyapunov(1)];
        fprintf('Counter: %.f and MLE: %f\n',counter,Lyapunov(1));
    end

    % Update counter
    counter = counter + 1;
end

%% Make the MLE vs ecc plot

% Alpha = 3/11 = 0.2727 exponents
alpha1_ecc = [0, 0.2096, 0.2917, 0.3513, 0.3997, 0.4401, 0.4750, 0.5056, 0.533];
alpha1_MLE = [0.0082413, 0.07916, 0.085404, 0.078165, 0.1211, 0.1191, 0.11701, 0.12315, 0.16284];

% Alpha = 3/11 = 0.1818 exponents
alpha2_ecc = [0, 0.2917, 0.3997, 0.4750, 0.533, 0.5797, 0.6186, 0.6515, 0.6799, 0.7045, 0.7262, 0.7361, 0.745];
alpha2_MLE = [0.0074608, 0.0058669, 0.0052082, 0.0060306, 0.0084947, 0.0070855, 0.0085782, 0.0091884, 0.01018, 0.01342, 0.013196, 0.019554, 0.09563];

% Alpha = 0.21 exponents
alpha3_ecc = [0, 0.2096, 0.2917, 0.3513, 0.3997, 0.4401, 0.4750, 0.5056, 0.5329, 0.5797, 0.6186, 0.6515, 0.6799];
alpha3_MLE = [0.0056548, 0.0042315, 0.0078868, 0.010472, 0.012174, 0.013114, 0.013777, 0.01423, 0.015689, 0.017453, 0.10612, 0.083915, 0.11612];

% Alpha = 0.24 exponents
alpha4_ecc = [0, 0.2096, 0.2917, 0.3513, 0.3997, 0.4401, 0.4750, 0.5056, 0.5329, 0.5797, 0.6186, 0.6515, 0.6799];
alpha4_MLE = [0.0099931, 0.015466, 0.020636, 0.097208, 0.10663, 0.10238, 0.090789, 0.083809, 0.093462, 0.12451, 0.13284];

figure (1)
plot(alpha1_ecc, alpha1_MLE,'-','LineWidth',2)
hold on
plot(alpha4_ecc, alpha4_MLE,'-','LineWidth',2)
plot(alpha3_ecc, alpha3_MLE,'-','LineWidth',2)
plot(alpha2_ecc, alpha2_MLE,'-','LineWidth',2)

plot(alpha1_ecc, alpha1_MLE,'.k','MarkerSize',10)
plot(alpha4_ecc, alpha4_MLE,'.k','MarkerSize',10)
plot(alpha3_ecc, alpha3_MLE,'.k','MarkerSize',10)
plot(alpha2_ecc, alpha2_MLE,'.k','MarkerSize',10)
hold off
grid on

xlabel('Eccentricity, $\varepsilon$','Interpreter','latex','FontSize', fs)
ylabel('MLE, $\lambda_1$','Interpreter','latex','FontSize', fs)
legend('$\alpha = 0.27$', '$\alpha = 0.21$', '$\alpha = 0.18$','Interpreter','latex',...
    'FontSize', fs,'location','northwest')

ax = gca;
ax.FontSize = fs-1;


%% Function definitions
% RHS of nonlinear isothermal coords ODE
function dydx = odefcn(theta, phi, a, R, r)
  gamma = sqrt((a+r*cos(theta)).^2.*sin(phi).^2 + (R+r*cos(theta)).^2.*cos(phi).^2);
  dydx = -1i*(r/gamma);
end
