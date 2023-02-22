% Numerical Jacobian matrix
% See (https://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf)
function J = myjacobian(func,x)
    % Set numerical derivative parameters
    N = length(x);

    %F_at_x = feval(func,x);
    epsilon = 1E-4;

    % Compute numerical derivative
    xperturb = x;
    xperturb2 = x; 
    xperturb_minus = x;
    xperturb_minus2 = x;
    %xperturb = x + epsilon;
    J = zeros(N);

    for i = 1:N
        xperturb(i) = xperturb(i) + epsilon; % x + h
        xperturb2(i) = xperturb2(i) + 2*epsilon; % x + 2h

        xperturb_minus(i) = xperturb_minus(i) - epsilon; % x - h
        xperturb_minus2(i) = xperturb_minus2(i) - 2*epsilon; % x - 2h
        
        % First order scheme
        %J(:,i) = (feval(func,xperturb) - F_at_x)/epsilon;

        % Second order scheme
        %J(:,i) = (feval(func,xperturb) - feval(func,xperturb_minus))/(2*epsilon);

        % Fourth order scheme
        J(:,i) = (-feval(func,xperturb2) + 8*feval(func,xperturb) ...
            - 8*feval(func,xperturb_minus) + feval(func,xperturb_minus2))/(12*epsilon);
        
        % Reset loop variables
        xperturb(i) = x(i);
        xperturb2(i) = x(i);
        xperturb_minus(i) = x(i);
        xperturb_minus2(i) = x(i);
    end
end
