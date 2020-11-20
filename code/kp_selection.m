function [x,y] = kp_selection(X,F)
%KP_SELECTION Selection operator
%
%   Inputs:
%   X - Solutions
%   F - Fitness values
%
%   Outputs:
%   x - First parent
%   y - Second parent

% Solution counter
k = size(X,1);

% Two random pairs
fp = randperm(k,2);
sp = randperm(k,2);

% New second pair if the former are equal
while sort(fp) == sort(sp)
    sp = randperm(k,2);
end

% Tournament selection
f1 = sum(F(fp'));
f2 = sum(F(sp'));

if f1 < f2
    x = X(fp(1),:);
    y = X(fp(2),:);
elseif f1 > f2
    x = X(sp(1),:);
    y = X(sp(2),:);
else
    r = rand;
    if r < 0.5
        x = X(fp(1),:);
        y = X(fp(2),:);
    else
        x = X(sp(1),:);
        y = X(sp(2),:);
    end
end

end