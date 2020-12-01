function [ND,D] = pareto_dominance(Y)
%PARETO_DOMINANCE Pareto dominance classification
%
%   Inputs:
%   Y - Objective values of solutions
%
%   Outputs:
%   ND - Non-dominated set
%   D - Dominated set

% Number of solutions
n = size(Y,1);

% Solutions
x = (1:n)';

% Dominated set
D = false(n,1);

% Main loop
for i = 1:n
    if D(i) == false
        for j = i+1:n
            if D(j) == false
                if D(j) == false && prod(Y(i,:)>=Y(j,:)) == 1 && sum(Y(i,:)>Y(j,:)) >= 1
                    D(j) = true;
                end
            end
            if prod(Y(j,:)>=Y(i,:)) == 1 && sum(Y(j,:)>Y(i,:)) >= 1
                D(i) = true;
                break;
            end
        end
    end
end

% Remove duplicates
D = find(D);

% Non-dominated set
ND = setdiff(x,D);

end