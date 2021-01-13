function y = ExplicitEuler(f, T, x0, N)
    % f: Given function
    % T: end time 
    % x0: initial value 
    % N: number of steps
    % Step 1: Discretize the domain
    t = linspace(0, T, N);
    % Step 2: Define time-step
    h = 1./N;
    % Step 3: Initialize the matrix holding the approximation of the
    % solution at each time step
    y = zeros(1, length(t));
    y(1) = x0;
    % Step 4: Implement the explicit euler method
    for n = 1:length(t)
        y(n+1) = y(n) + h.*f(y(n));
    end