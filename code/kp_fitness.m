function F = kp_fitness(Z)
%KP_FITNESS Calculate fitness values
%
%   Inputs:
%   Z - Solutions
%
%   Outputs:
%   F - Fitness assignment

% Dimensions
[n,p] = size(Z);

% Fitness
F = zeros(n,1);

% Solutions to rank
Y = Z;

% Domination counter
dc = 1;

% Ranked counter
rc = 0;

% First front
[ND,~] = pareto_dominance(Y);
F(ND) = dc;
Y(ND,:) = -Inf*ones(length(ND),p);
dc = dc + 1;
rc = rc + length(ND);

% Following fronts
while true
    if rc == n
        break;
    end
    [ND,~] = pareto_dominance(Y);
    F(ND) = dc;
    Y(ND,:) = -Inf*ones(length(ND),p);
    dc = dc + 1;
    rc = rc + length(ND);
end

end