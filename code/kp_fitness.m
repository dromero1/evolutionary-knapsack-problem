function F = kp_fitness(Z)
%KP_FITNESS Calculate fitness values
%
%   Inputs:
%   Z - Solutions
%
%   Outputs:
%   F - Fitness assignment

% Dimensions
[k,p] = size(Z);

% Fitness
F = zeros(k,1);

% Solutions to rank
Y = Z;

% Domination counter
dc = 1;

% Ranked counter
rc = 0;

% Unfeasible solutions
Iuf = find(Y(:,end)~=1);
Y(Iuf,:) = (-1)*ones(length(Iuf),p);

% First front
[ND,~] = pareto_dominance(Y);
F(ND) = dc;
Y(ND,:) = (-Inf)*ones(length(ND),p);
dc = dc + 1;
rc = rc + length(ND);

% Following fronts
while rc ~= k
    [ND,~] = pareto_dominance(Y);
    F(ND) = dc;
    Y(ND,:) = (-Inf)*ones(length(ND),p);
    dc = dc + 1;
    rc = rc + length(ND);
end

end