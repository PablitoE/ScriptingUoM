% Autofocusing functions for blockproc (except CC)
function funs = autofocusing_funs()
funs.Entropy = @(block) entropy(block.data);
funs.L1 = @(block) mean(block.data(:));
funs.Gradient = @(block) get_gradient(block.data);
funs.Laplacian = @(block) get_laplacian(block.data);
funs.Variance = @(block) var(block.data(:));
funs.Tamura = @(block) sqrt(std(block.data(:))/mean(block.data(:)));
funs.CC = @CC_in;
funs.RC = @(block) get_RC(block.data);
funs.Gini = @(block) get_Gini(block.data);
funs.StdGrad = @(block) get_StdGrad(block.data);
end

function gra = get_gradient(A)
% Gradient
GRA = cat(3, A(2:end, 2:end) - A(1:end-1, 2:end), A(2:end, 2:end) - A(2:end, 1:end-1));
gra = mean2(sqrt(sum(GRA.^2, 3)));
end

function lap = get_laplacian(A)
% Laplacian
LAP = A(1:end-2, 2:end-1) + A(3:end, 2:end-1) + A(2:end-1, 1:end-2) + A(2:end-1, 3:end) - 4 * A(2:end-1, 2:end-1);
lap = mean2(LAP.^2);
end

function Bcc = CC_in(A,B,w)
[M, N] = size(A);
Mo = floor(M/w);
No = floor(N/w);
Bcc = zeros(Mo,No);
for kf = 1:Mo
   for kc = 1:No
      vA = A((kf-1)*w+1:kf*w,(kc-1)*w+1:kc*w);
      vB = B((kf-1)*w+1:kf*w,(kc-1)*w+1:kc*w);
      R = corrcoef(vA(:), vB(:));
      Bcc(kf,kc) = R(2);
   end
end
end
      
function rc = get_RC(Psi)
% RC method
var_real = var(real(Psi(:)));
var_imag = var(imag(Psi(:)));
rc = var_real / var_imag;
end

function gini = get_Gini(A)
% Gini
N = numel(A);
A = sort(A(:));
L1 = mean(A(:));
gini = 1 - sum((N - (1:N)' + 0.5) .* A) * 2 / L1 / N^2;
end

function std = get_StdGrad(A)
[gx, gy] = gradient(A);
grad = hypot(gx, gy);
grad = grad(~isnan(grad));
meany = mean(grad);
std = sqrt(sum((grad - meany).^2));
end