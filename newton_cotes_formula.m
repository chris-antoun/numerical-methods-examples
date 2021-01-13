% This script implements the Newton-Cotes formula
% Define a function called newton_cotes_formula
function y = newton_cotes_formula(f,x_i,x)
    % f is the function we are trying to approximate
    % x_i is a vector containing the interpolating points
    % x is a vector of points that we want to interpolate over
    a = x(:,1);
    b = x(:,end);
    L = @(x) 1;
    y = 0;
    for j = 1:length(x_i)
        for i = 1:length(x_i)
            if i ~= j
                L = @(x) L(x).*(x-x_i(j))./(x_i(i)-x_i(j));
            end
        end
    end
    
    for j = 1:length(x)
        y = y+ f(x(j)).*integral(L, a, b); 
    end

end
