function zk = newtonODE(f, df, alpha ,beta, N, tolerance)
    % f: given function; input to compute the vector b
    % df: function; the derivative of the function f;
    % alpha: scalar; input to compute the vector b 
    % beta: scalar; input to compute the vector b
    % N: scalar; input to determine the dimension 
    % tolerance: scalar; input to determine when the approximated solution
    % is "close" enough
    % Step 0: Compute the matrix A
    % Step 0.1. Initialize A
    A = zeros(N+1, N+1);
    % Step 0.2. Define the first row of A
    A(1,1) = 1;
    % Step 0.3. Define the last row of A
    A(N+1, N+1) = 1;
    % Step 0.4. Define A
    for row = 2:N
        for col = 0:row-2
            row_to_add = [zeros(1, col) 1 -2 1 zeros(1, N+1-3-col)];
            A(row, :) = row_to_add;
        end
    end
    % Step 1: Initialize the matrix z
    z0 = zeros(N+1,1);
    z = [z0, zeros(N+1,1)];
    % Step 2: Initialize the error and the running index 
    error = tolerance + 1;
    k = 0;
    % Step 3: Compute zk
    while error > tolerance
        k = k + 1;
        % Step 3.1: Compute zk
            % Step 3.1.1 Define b
            b = zeros(N+1,1);
            % Step 3.1.2. Define first entry of b
            b(1) = alpha;
            % Step 3.1.3. Define last entry of b
            b(N+1) = beta;
            % Step 3.1.4. Define b
            for row = 2:N
                b(row) = (1./N).^2.*f(z(row,k));
            end
        % Step 3.2 Define g
        g = A*z(:,k) - b;
        % Step 3.3. Compute the gradient of g
        Dg = A; % Initializing Dg to be A as it is the same except on the diagonal
        for i = 2:N
            Dg(i,i) = Dg(i,i) - df(z(i,k));
        end
        % Step 3.4. Compute zk
        z(:,k+1) = z(:,k) - Dg\g;
        % Step 4: Compute the error
        % Step 4.1. Compute the sum
        summand = 0;
        for i = 1:N+1
            summand = summand + abs(z(i,k+1) - z(i,k)).^2;
        end
        % Step 4.2. Compute the error
        error = sqrt(1./N.*summand);
    end       
    zk = z(:, end);
    