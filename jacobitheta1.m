function z = jacobitheta1(w,p,cap)
% This function computes the first Jacobi theta function, where:
% z = complex number output
% w = u + iv = complex number input
% p = e^(i*pi*t) is the nome
% cap = truncation error of the infinite sum

u = real(w);
v = imag(w);

sum = 0;
for n = 0:cap
    sum = sum + (-1)^n*p^(n*(n+1))*sin((2*n+1).*(u+1i.*v));
end

z = 2*p^(1/4)*sum;

end