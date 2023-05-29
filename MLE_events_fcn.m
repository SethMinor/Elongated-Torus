function [position,isterminal,direction] = MLE_events_fcn(~,y,delta,ExpFac)
  y1 = y(1:4);
  y2 = y(5:8);
  position = norm(y1-y2) - delta*ExpFac; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 1;   % Stop when norm is increasing
end
