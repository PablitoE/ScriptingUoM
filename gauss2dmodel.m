function y = gauss2dmodel(x,xdata)
   m = x(1:2);
   cov = [[x(3), x(4)];[x(4),x(5)]];
   a = x(6);
   xdata = xdata - m;
   y = a*exp(-0.5*sum((xdata/cov).*xdata,2));
end