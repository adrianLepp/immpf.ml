function [x_k,y_k] = solveThreeTank(x_k,parameter)
%solveDreitank: solve equations of three tank system in discrete form for
%one time step by approximation of Runge Kutta method 4. order

%___INPUTS:
%   x_k         state x for k-1
%   parameter   parameters for the system equations, noise sigma , time step dt

%___OUTPUTS:
%   x_k         state x for k

        u = parameter.u;
        c13 = parameter.c13;
        c32 = parameter.c32;
        cA2 = parameter.cA2;
        A = parameter.A;
        g = parameter.g;

        c = parameter.c; %measurement matrix: y=h(x)=c*x
        
        %Noise
        sigmaX = parameter.sigmaX; % process noise
        sigmaY = parameter.sigmaY; % measurement noise
     
        T = parameter.dt; % Time Step k
        dt = T/2; % solve in 2 steps

        % Runge-Kutta 4. order
        for tau = dt : dt : T
            xtemp = x_k;
            dx1(1,1) = 1/A*(u-c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3))));
            dx1(2,1) = 1/A*(c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2)))-cA2*sqrt(2*g*abs(xtemp(2))));
            dx1(3,1) = 1/A*(c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3)))-c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2))));
            dx1 = dx1 * dt;
            xtemp = x_k + dx1 / 2;
            dx2(1,1) = 1/A*(u-c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3))));
            dx2(2,1) = 1/A*(c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2)))-cA2*sqrt(2*g*abs(xtemp(2))));
            dx2(3,1) = 1/A*(c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3)))-c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2))));
            dx2 = dx2 * dt;
            xtemp = x_k + dx2 / 2;
            dx3(1,1) = 1/A*(u-c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3))));
            dx3(2,1) = 1/A*(c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2)))-cA2*sqrt(2*g*abs(xtemp(2))));
            dx3(3,1) = 1/A*(c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3)))-c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2))));
            dx3 = dx3 * dt;
            xtemp = x_k + dx3;
            dx4(1,1) = 1/A*(u-c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3))));
            dx4(2,1) = 1/A*(c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2)))-cA2*sqrt(2*g*abs(xtemp(2))));
            dx4(3,1) = 1/A*(c13*sign(xtemp(1)-xtemp(3))*sqrt(2*g*abs(xtemp(1)-xtemp(3)))-c32*sign(xtemp(3)-xtemp(2))*sqrt(2*g*abs(xtemp(3)-xtemp(2))));
            dx4 = dx4 * dt;

            x_k = x_k + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6 + randn(3,1)*sqrt(sigmaX);

            y_k = c * x_k + randn(3,1)*sqrt(sigmaY);
        end
end
