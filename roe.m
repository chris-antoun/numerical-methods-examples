% clear workspace 
clear
clc

% define variables
xmin = -2;
xmax = 3; 
tmin = 0;
tmax = 1;

% discretize the domain
dx = 1/32;
x = xmin:dx:xmax;
dt = 1/2*dx*1/2;

% set initial condition
u0(x < 0) = -2;
u0(x > 0 & x < 1) = 2;
u0(x > 1) = 0;
u = u0;
unp1 = u;

% define flux
f = @(u) 1/2.*u.^2;

% define numerical flux at time = 0
N = (xmax - xmin)/dx;
df = @(u) u;
roe_average = zeros(size(x));
F = zeros(size(x));

% calculate the initial roe average

 for i = 1:N
     if u(i+1) ~= u(i)
        roe_average(i) = (f(u(i+1)) - f(u(i)))/(u(i+1)-u(i));
    else
        roe_average(i) = df(u(i));
    end
 end
 
 % calculate the initial flux 
 
 for i = 1:N-1
    if roe_average(i) >= 0
        F(i) = f(u(i));
    else
        F(i) = f(u(i+1));
    end
end

% loop through time
t = tmin;
nsteps = tmax/dt;

for n = 1 : nsteps
    
    for i = 1 : N
       
        if x(i) == xmin || x(i) == xmax
            % implement neumann boundary condition
            unp1 = u;
        else
            % implement roe scheme
            if u(i+1) ~= u(i)
                roe_average(i) = (f(u(i+1)) - f(u(i)))/(u(i+1)-u(i));
            else
                roe_average(i) = df(u(i));
            end
        
        if roe_average(i) >= 0
            F(i) = f(u(i));
        else
            F(i) = f(u(i+1));
        end
            unp1(i) = u(i) - dt/dx.*(F(i)-F(i-1));
        end
 
    end
    % calculate exact solution
    exact = entropy_solution(x,t);
     
    % update t and u
    t = t + dt;
    u = unp1;    
    
    % plot the solution
    plot(x, u, 'b-');
    hold on
    plot(x, exact, 'r-')
    hold off
    xlim([xmin xmax])
    legend('numerical','exact')
    xlabel('x', 'fontsize', 16)
    ylabel('U(t,x)','fontsize',16)
    title(sprintf('time = %1.3f',t), 'fontsize',16)
    grid
    shg
    pause(dt)
    
end


function exact = entropy_solution(x,t)
    for i = 1:length(x)
        if x(i) <= -2*t
            exact(i) = -2;
        else 
            if x(i) <= 2*t
                exact(i) = x(i)/t;
            else 
                if x(i) <= 1 + t
                    exact(i) = 2;
                else 
                    exact(i) = 0;
                end
            end
        end
    end
end
    