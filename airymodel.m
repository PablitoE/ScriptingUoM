% Airy amplitude : a*J1(sqrt(x'cov-1x)) / x'cov-1x, J1 is the Bessel function of the first kind of order one
function y = airymodel(x,xdata)
   m = x(1:2);
   cov = [[x(3), x(4)];[x(4),x(5)]];
   a = x(6);
   xdata = xdata - m;
   x_eval = sqrt(sum((xdata/cov).*xdata,2));
   y = zeros(size(x_eval));
   x_zero = x_eval == 0;
   y(x_zero) = a;
   x_eval = x_eval(~x_zero);
   y(~x_zero) = 2*a/1.22*besselj(1,x_eval*1.22)./x_eval;
end