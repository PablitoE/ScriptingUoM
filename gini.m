function gini_coef=gini(values, pop)
% FORMAT: gini(values, population)
% Givne these values, this computes the GINI coefficient according to
% footnote 17 of Noorbakhsh
if size(values,2)~=1 || size(pop,2)~=1 || size(values,1)~=size(pop,1)
    error('Arguments need to be column vectors')
end
P = sum(pop);
mu = sum(values.*(pop/P));
C = size(values,1);
total = 0;

for i=1:C
    for j=1:C
        total=total+pop(i)*pop(j)/(P^2)*abs(values(i)-values(j));
    end
end

gini_coef = total / mu;