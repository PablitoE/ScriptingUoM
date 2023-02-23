function root = get_ASM_root(dx,dy,Nc,Nr,wavel)

dx = single(dx); dy = single(dy);
Nc = single(Nc); Nr = single(Nr);
wavel = single(wavel);

% Find the point of nux and nuy at which the propagator becomes
% undersampled at distance maxz. A little work will show this to be when
% nuy = 0, nux = +/- (lambda * sqrt((2 maxz dnux)^2 + 1 ))^-1
dnux  = 1/(dx * Nc);   % Frequency 'pixel width'
dnuy  = 1/(dy * Nr);

x = - Nc/2:Nc/2-1;  % xs and ys
y = - Nr/2:Nr/2-1;
[xx,yy] = meshgrid(x,y);        % root in phase multiplier
clear('x', 'y');

root = wavel^2 *( (xx.*dnux).^2 + (yy.*dnuy).^2);
root(root > 1) = 1;
root = abs(sqrt(1 - root));