function [majsiz, minsiz] = gaussfit_im(im)
   [m,n] = size(im);
   % Remove nans from the analysis
   [valid_rows, valid_cols] = find(~isnan(im));
   valid_lin= sub2ind([m,n],valid_rows,valid_cols);
   % Prepare inputs to lsqnonlin
   y = im(valid_lin);   % data
   [cc, rr] = meshgrid(1:n, 1:m);
   x0 = [m/2,n/2,1,0,1,max(y)]; % [mean_col, mean_row, var_col, cov, var_row, amplitude] 
   xdata = [cc(valid_lin), rr(valid_lin)];   % col is xx, row is yy
   % Bounds
   lb = [0, 0, 0, -m*n, 0, -inf];
   ub = [n+1, m+1, n^2, m*n, m*2, inf];
   options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','none');
   x = lsqcurvefit(@gauss2dmodel,x0,xdata,y,lb,ub,options);
   cov = [[x(3), x(4)];[x(4),x(5)]];
   e = eig(cov);
   majsiz = e(2);
   minsiz = e(1);
end