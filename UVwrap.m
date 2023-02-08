% Wrap U,V to interval
function wrapped = UVwrap(array, interval)
    wrapped = mod(array - interval(1), range(interval)) + interval(1);
end
