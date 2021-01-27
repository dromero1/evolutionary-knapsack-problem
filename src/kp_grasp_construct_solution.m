function [x,fea,iter] = kp_grasp_construct_solution(n,m,W,A,b,alpha)
%KP_GRASP_CONSTRUCT_SOLUTION Constructive phase of KP GRASP approximation
%
%   Inputs:
%   n - Number of items
%   m - Number of constraints
%   W - Objective coefficients
%   A - Constraint coefficients
%   b - Resource capacity
%   alpha - Best candidate percentage
%
%   Outputs:
%   x - Solution
%   fea - Feasibility percentage
%   iter - Iterations

%% Initialization phase

% 1 if element i is added to the knapsack
x = true(n,1);

% Resource consumption
R = A*x;

% Maximum number of iterations
max_iter = 10^3;

%% Feasibility lookup phase

% Positive constraints
Ip = (b>=0);

% Negative constraints
In = (b<0);

% Inverse of mean benefits
mW = 1./mean(W);

% Selection criteria for removal
if sum(Ip) > 1
    mur = sum(A(Ip,:))';
else
    mur = A(Ip,:)';
end

% Selection criteria for addition
mua = mean(A(In,:))';

% Iterations
iter = 0;

while prod(R <= b) ~= 1 && iter < max_iter
    if any(x)
        % Selection criteria for removal
        mu = mur.*x;
        % Selection threshold
        if alpha ~= 0
            mimu = min(mu(x));
            th = (1-alpha)*(max(mu(x))-mimu)+mimu;
        else
            th = max(mu);
        end
        RCL = (mu>=th);
        % Select element to remove
        if sum(RCL) > 1
            elm = find(RCL);
            p = mW(RCL)./sum(mW(RCL));
            [sP,Isp] = sort(p,'descend');
            r = rand;
            cP = 0;
            for i = 1:length(elm)
                cP = cP + sP(i);
                if r <= cP
                    idx = Isp(i);
                    k = elm(idx);
                    break;
                end
            end
        else
            k = find(RCL);
        end
        % Remove selected element from knapsack
        x(k) = 0;
        % Update resource consumption
        R = R - A(:,k);
    end
    if prod(R(In) <= b(In)) ~= 1
        % Selection criteria for addition
        y = ~x;
        mu = mua.*(y);
        % Selection threshold
        mimu = min(mu(y));
        th = alpha*(max(mu(y))-mimu)+mimu;
        RCL = (mu<=th);
        % Select element to add
        elm = find(RCL);
        idx = randi(length(elm),1);
        k = elm(idx);
        % Add selected element to knapsack
        x(k) = 1;
        % Update resource consumption
        R = R + A(:,k);
    end
    iter = iter + 1;
end

% Feasibility percentage
fea = sum(A*x <= b)/m;

end