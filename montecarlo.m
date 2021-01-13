function appr = montecarlo(f,M)
    % f: the function we want to approximate its integral
    % M: the number of uniformly distributed r.v.s chosen to be
    % interpolating points
    
    % Step 1: Pick M random numbers in the domain
    interpolating_points = rand(M,1);
    % Step 2: Define the approximated integral
    % Initialize approximated integral
    appr = 0;
    for m = 1:M
        appr =  appr +(1./M).*f(interpolating_points(m));
    end