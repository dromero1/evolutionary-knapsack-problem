function [X,Z,nsol] = kp_ga(ti,n,p,m,W,A,b,alpha,J,mup,mt,dbg)
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
%   mup - Mutation probability
%   mt - Maximum execution time
%   dbg - Debug mode
%
%   Outputs:
%   X - Solutions
%   Z - Objective values
%   nsol - Number of solutions

% Initial time
t0 = toc;

% Number of generations
num_gen = 100;

% Population size
pop_size = 500;

% Number of children
num_children = 100;

% Solutions
X = false(pop_size,n);
Z = zeros(pop_size,p+1);

% Number of feasible solutions
fc = 0;

% Number of solutions
nsol = 0;

% Generate initial population
for i = 1:pop_size
    [x,fea,~] = kp_grasp_construct_solution(n,m,W,A,b,alpha);
    X(i,:) = x';
    Z(i,:) = [(W*x)' fea];
    nsol = nsol + 1;
end
% Initial fitness assignment
F = kp_fitness(Z);
% Evolution
for gen = 1:num_gen
    % Check time
    if toc - t0 > mt
        break;
    end
    % Children solutions
    S = false(num_children,n);
    Zs = zeros(num_children,p+1);
    % Genetic operators
    for j = 1:num_children
        % Check time
        if toc - t0 > mt
            break;
        end
        % Selection
        [p1,p2] = kp_selection(X,F);
        % Crossover
        ch = kp_crossover(p1,p2);
        % Feasibility percentage
        ch_fea = sum(A*ch' <= b)/m;
        % Local search
        if ch_fea == 1
            ch_star = kp_vnd(ch',n,m,W,A,b,J,false,t0,mt);
            ch = ch_star;
        end
        % Mutation
        r = rand;
        if r < mup
            ch = kp_mutation(ch,mup);
            ch_fea = sum(A*ch' <= b)/m;
        end
        % Save child
        S(j,:) = ch;
        % Save objetive values
        Zs(j,:) = [(W*ch')' ch_fea];
        % Update solution count
        nsol = nsol + 1;
    end
    % Merge solutions
    X = [X; S];
    Z = [Z; Zs];
    % Fitness assignment
    F = kp_fitness(Z);
    % Update population
    [~,Ipu] = sort(F);
    Ipu = Ipu(1:pop_size);
    X = X(Ipu,:);
    Z = Z(Ipu,:);
    F = F(Ipu,:);
    if dbg == true
        fprintf('GA Instance %d (gen. = %d, mup. = %0.2f, ',ti,gen,mup);
        fprintf('fitness std. = %0.2f)\n',std(F));
    end
end
% Improvement phase
for i = 1:pop_size
    % Check time
    if toc - t0 > mt
        break;
    end
    % Get solution
    x = X(i,:);
    % Get feasibility percentage
    fea = Z(i,p+1);
    if fea == 1
        % Update feasible count
        fc = fc + 1;
        % Variable neighborhood descent
        x_star = kp_vnd(x',n,m,W,A,b,J,false,t0,mt);
        % Update solution
        X(i,:) = x_star;
        % Update objetive values
        Z(i,:) = [(W*x_star')' fea];
    end
end

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