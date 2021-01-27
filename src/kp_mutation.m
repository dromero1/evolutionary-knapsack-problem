function x = kp_mutation(x,gamma)
%KP_MUTATION Mutation operator
%
%   Inputs:
%   x - Solution
%   gamma - Bitwise mutation probability
%
%   Outputs:
%   x - Mutated solution
%
%   Idea taken from:
%   Eiben, A. E., & Smith, J. E. (2015). Introduction to evolutionary
%   computing. Springer. (Pag. 52)

% Items count
n = length(x);

% Bitwise mutation
for i = 1:n
    r = rand;
    if r < gamma
        x(i) = ~x(i);
    end
end

end