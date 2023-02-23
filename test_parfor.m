% Trying a progress bar in a parfor loop
% if isempty(gcp("nocreate"))
%    parpool(3)
% end
fprintf('Progress:\n');
m = 10;
fprintf(['\n' repmat('.',1,m) '\n\n']);
parfor k = 1:m
   time_to_suspend = 1 + rand() * 2;
   pause(time_to_suspend)
   fprintf('|\n');
end
% delete(gcp('nocreate'));