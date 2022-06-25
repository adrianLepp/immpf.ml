function [xK,yK] = solveThreeTank(xK,parameter)
%solveDreitank: solve equations of three tank system in discrete form for
%one time step by approximation of Runge Kutta method 4. order

%25.06.2022 Adrian Lepp

%___INPUTS:
%   xK          state x for k-1
%   parameter   parameters for the system equations, noise sigma , time step dt

%___OUTPUTS:
%   xK          state x for k
%   yK          measurement y for k
 
        c = parameter.c; % measurement matrix: y=h(x)=c*x
        
        %Noise
        sigmaX = parameter.sigmaX; % process noise
        sigmaY = parameter.sigmaY; % measurement noise
     
        T = parameter.dt; % Time Step k
        dt = T/2; % solve in 2 steps

        % Runge-Kutta 4. order
        for tau = dt : dt : T
            xtemp = xK;
            dx1 = solveThreeTankLinear(xtemp, parameter)*dt;
            xtemp = xK + dx1 / 2;
            dx2 = solveThreeTankLinear(xtemp, parameter)*dt;
            xtemp = xK + dx2 / 2;
            dx3 = solveThreeTankLinear(xtemp, parameter)*dt;
            xtemp = xK + dx3 / 2;
            dx4 = solveThreeTankLinear(xtemp, parameter) *dt;

            xK = xK + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
        end
        xK = xK + randn(3,1) * sqrt(sigmaX);
        yK = c * xK + randn(3,1) * sqrt(sigmaY);
        
    function [dx] =  solveThreeTankLinear(x, parameter)
            u = parameter.u;
            c13 = parameter.c13;
            c32 = parameter.c32;
            cA2 = parameter.cA2;
            A = parameter.A;
            g = parameter.g;
            dx = zeros(3,1);
            
            dx(1) = 1/A*(u-c13*sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3))));
            dx(2) = 1/A*(c32*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2)))-cA2*sqrt(2*g*abs(x(2))));
            dx(3) = 1/A*(c13*sign(x(1)-x(3))*sqrt(2*g*abs(x(1)-x(3)))-c32*sign(x(3)-x(2))*sqrt(2*g*abs(x(3)-x(2))));
    end
end
