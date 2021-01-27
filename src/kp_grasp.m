function [X,Z,nsol] = kp_grasp(ti,n,p,m,W,A,b,alpha,J,mt,dbg,ls)
%KP_GRASP GRASP method approximation to the knapsack problem
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
%   ls - Local search mode
%
%   Outputs:
%   X - Solutions
%   Z - Objective values
%   nsol - Number of solutions

% Initial time
t0 = toc;

% Solutions
X = false(1,n);
Z = zeros(1,p+1);

% Number of feasible solutions
fc = 0;

% Main loop
i = 1;
while toc - t0 <= mt
    % Randomized constructive solution
    [x,fea,~] = kp_grasp_construct_solution(n,m,W,A,b,alpha);
    if ls == true
        if fea == 1
            % Variable neighborhood descent
            X_star = kp_vnd(x,n,m,W,A,b,J,true,t0,mt);
            % Save local search solutions
            n_star = size(X_star,1);
            for j = 1:n_star
                % Local search solution
                x_star = X_star(j,:);
                X(i,:) = x_star;
                % Determine feasibility
                fea = sum(A*x_star' <= b)/m;
                % Local search objective values
                Z(i,:) = [(W*x_star')' fea];
                % Update feasible count
                if fea == 1
                    fc = fc + 1;
                end
                % Display
                if dbg == true
                    fprintf('GRASP Instance %d (alpha = %0.2f, ',ti,alpha);
                    fprintf('rep. = %d, feas. = %0.2f)\n',i,fea);
                end
                i = i + 1;
            end
        else
            % Save unfeasible solution
            X(i,:) = x;
            Z(i,:) = [(W*x)' fea];
            i = i + 1;
        end
    else
        % Save solution
        X(i,:) = x;
        Z(i,:) = [(W*x)' fea];
        if fea == 1
            fc = fc + 1;
        end
        % Display
        if dbg == true
            fprintf('GRASP Instance %d (alpha = %0.2f, ',ti,alpha);
            fprintf('rep. = %d, feas. = %0.2f)\n',i,fea);
        end
        i = i + 1;
    end
end

% Number of solutions
nsol = i - 1;

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