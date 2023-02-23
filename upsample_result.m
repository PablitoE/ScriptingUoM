% Upsample by a factor of x means getting x samples for each input sample
function fout = upsample_result(f,factor,desired_size)
[Nr, Nc] = deal(desired_size(1), desired_size(2));
F = griddedInterpolant(f);
cc = linspace(1,(Nc-1)/factor+1,Nc)';
rr = linspace(1,(Nr-1)/factor+1,Nr)';
fout = F({rr,cc});
