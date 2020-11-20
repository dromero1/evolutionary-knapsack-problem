function [X,Z,nsol] = kp_ga(ti,n,p,m,W,A,b,alpha,J,mt,dbg)
%KP_GA Genetic approach to the knapsack problem
%
%   Inputs:
%   ti - Test instance
%   n - Number of items
%   p - Number of objectives
%   m - Number of constraints
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   alpha - Best candidate percentage
%   J - Number of neighborhoods
%   mt - Maximum execution time
%   dbg - Debug mode
%
%   Outputs:
%   X - Solutions
%   Z - Objective values
%   nsol - Number of solutions

% Initial time
t0 = toc;

% Population size
pop_size = 100;

% Number of children
num_children = 30;

% Mutation probability
epsilon = 0.1;

% Solutions
X = false(pop_size,n);
Z = zeros(pop_size,p+1);
fc = 0;

% Generate initial population
for i = 1:pop_size
    [x,fea,~] = kp_grasp_construct_solution(n,m,W,A,b,alpha);
    X(i,:) = x';
    Z(i,:) = [(W*x)' fea];
end

% Main loop
while toc - t0 <= mt
    % Fitness assignment
    F = kp_fitness(Z);
    % Genetic operators
    for j = 1:num_children
        
    end
    % Update population
end

% Variable neighborhood search

% Remove duplicates
[X,ix,~] = unique(X,'rows');
Z = Z(ix,:);

% Remove infeasible solutions
if fc >= 1
    If = (Z(:,p+1)==1);
    X = X(If,:);
    Z = Z(If,:);
end

% Get non-dominated solutions
[Ipo,~] = pareto_dominance(Z);
X = X(Ipo,:);
Z = Z(Ipo,:);

end