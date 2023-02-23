syms fx fy lam d fm c dfx z

hz = 2*pi*(sqrt(lam^(-2) - fx^2 - fy^2) - c * fm - d/fm*(fx^2 + fy^2)) * z;
% The phase difference between samples is modelled with dhzdfx * dfx
dhzdfx = diff(hz, fx);
% Looking for a maximum
d2hzdfx2 = diff(dhzdfx, fx);
max_dhzdfx_expression = d2hzdfx2 == 0;
fx_max_dhzdfx_array = solve(max_dhzdfx_expression, fx);
fx_max_dhzdfx = fx_max_dhzdfx_array(6);   % Keep the real one
% Find the phase difference at the maximum
dhzdfx_max = subs(dhzdfx, fx, fx_max_dhzdfx);
% Simplifying fy = 0
simple_fx_max_dhzdfx = subs(fx_max_dhzdfx, fy, 0);
simple_dhzdfx_max = subs(dhzdfx_max, fy, 0);
% Well sampled hz condition
well_sampled_condition = abs(simple_dhzdfx_max) * dfx / pi == 1;
% Test using kmax = km == km di / (4 d z) sqrt(Nc"^2 + Nr"^2) : known d
using_fm = solve(well_sampled_condition, fm);