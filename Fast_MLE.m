%% Computes fast MLE
% Follows an IC and a nearby perturbation, iteratively rescaling the
% perturbation.

% NOTE FOR USE
% 'Dginv' is required, so isothermal coords must exist in workspace before 
% calling this function.

% INPUTS
% delta = perturbation distance (1E-6)
% ExpFac = inter-orbit distances >= ExpFac*delta are rescaled (1E2)
% y0 = [u1,u2,v1,v2] (row vector), initial condtion of interest
% tf = final integration time
% display = '1' (show plots) or '0' (no plots)
% params = {N,q,p,c,r,a,R,cap,phi,theta,gr,Dginv} (cell type)

% OUTPUTS
% MLE = largest LE, computed via median filtering with outlier removal
% Associated plots (for 'display = 1')

function [MLE] = Fast_MLE(delta, ExpFac, y0, tf, display, params)
% Initial random perturbation
rand_dir = rand(4,1); % Random direction
pert = delta*rand_dir/norm(rand_dir); % Perturbation
y1_0 = y0'; % First IC
y2_0 = y1_0 + pert; % Second (perturbed) IC
y0 = [y1_0; y2_0];
ye = y0;

% Main loop
t_cumulative = 0;
counter = 1;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'Events',...
    @(t,y) MLE_events_fcn(t,y,delta,ExpFac));

while t_cumulative < tf
    % Update IC with rescaled perturbation
    y0(1:4) = ye(1:4);
    pert_dir = ye(5:8) - ye(1:4); % Rescale in direction of 2nd orbit
    pert = delta*pert_dir/norm(pert_dir); % Perturbation
    y0(5:8) = y0(1:4) + pert;

    % Numerical integration using ode45
    [~,~,te,ye,~] = ode45(@(t,y) two_IC_torus(t,y,params),...
        [0, tf], y0, options);
    ye = ye';
    
    % Check for events
    if isempty(te) == 0
        % Add MLE to list
        lambdas(counter) = (1/te)*log(ExpFac);
    
        % Update loop variables
        t_cumulative = t_cumulative + te;
        times(counter) = t_cumulative;
        counter = counter + 1;
    end
end

% Average the MLEs (with filtering for spikes)
MLE = mean(medfilt1(rmoutliers(lambdas)));

% Optionally, display plots
if display == 1
    fs = 18; % Font size, for plotting

    figure (98)
    sgtitle("MLE $=$ "+MLE,'interpreter','latex','fontsize',fs)
    subplot(2,1,1)
    plot(times,lambdas)
    hold on
    plot(times,medfilt1(lambdas))
    hold off
    xlabel('Time, $t$','interpreter','latex','fontsize',fs)
    ylabel('Approx. $\lambda$','interpreter','latex','fontsize',fs)
    yline(MLE,'--r')
    legend('Unfiltered','Filtered','Computed',...
        'interpreter','latex','fontsize',fs)
    
    subplot(2,1,2)
    histogram(lambdas)
    hold on
    xline(MLE,'--r')
    xlabel('Approx. $\lambda$','interpreter','latex','fontsize',fs)
    ylabel('Counts','interpreter','latex','fontsize',fs)
    hold off
    legend('','Computed','interpreter','latex','fontsize',fs)
end
end
