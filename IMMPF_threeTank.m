% hybrid state estimation with IMMPF based on Blom
% realization for a simulated tree tank system

% Author: Adrian Lepp
% last change: 25.02.2022

% makes use of immpfCycle.m, solveThreeTank.m and dreitTank.mat 


%% load parameters
close all
clear
clc

load 'dreiTank.mat' %system parameters
n = parameter.n; %Dimension of x
dt = parameter.dt; %time step of k

%% init

x0 = zeros(n,1);    % initial state values     !!!
q = 4;              % Modes                     !!!
S = 100;            % Particles per Mode        !!!
T = 150;            % Simulation Time (s)       !!!
t = T/dt;           % time steps 

% Initialize modespecific values for system parameters !!!
parameterM(1:q) = struct(parameter);
parameterM(2).u = 1.5*parameter.u;
parameterM(3).c13 = 0.5*parameter.c13;
parameterM(4).c32 = 0.7*parameter.c32;

% Transition Probability Matrix
Pi = ones(q,q);
trans1 = 0.991; % i=j                       !!!
trans2 = (1 - trans1)/(q-1); % i!=j
for i = 1 : q
    for j = 1 : q
        if j == i
            Pi(i,j) = trans1;
        else
            Pi(i,j) = trans2;
        end
    end
end

% initial posterior mode probability
mInit(1)=trans1;            %!!!
mInit(2:q)=trans2;          %!!!

xPost = x0 .* ones(n,S,q) + sqrt(parameter.sigmaX) * randn(n,S,q);  % initial state particles
wPost = ones(S,q);                                                  % initial particle weights
for j = 1 : q
    wPost(:,j)  = wPost(:,j) .* mInit(j)./S;
end
% Output values that are stored over the whole time
x = zeros(n,t);         % real state
xEstM = zeros(n,q,t);   % mode specific state estimation
xEst = zeros(n,t);      % state estimation
mPost = zeros(q,t);     % posterior mode probability
    
%% Simulation
for k = 1 : t
    
    %system behaviour !!! -> time in seconds/dt
    if k <= 500
        %default system behaviour;
    elseif k > 500 && k <= 1000
        %fault in control input u equal to u(m2)
        parameter = parameterM(2);
    elseif k > 1000 && k <= 1500
        %default system behaviour;
        parameter = parameterM(1);
    end

    % solve equations f() and h() for real system
    if k == 1
        [x(:,k),y] = solveThreeTank(x0,parameter);
    else
        [x(:,k),y] = solveThreeTank(x(:,k-1),parameter);
    end

    % IMMPF for all q Modes and S Particles 
    [xPost,wPost] = immpfCycle(Pi,xPost,wPost,y,@solveThreeTank, parameterM,q,S,n);
    
    % Output
    for j = 1 : q
        % posterior mode probability
        mPost(j,k) = sum(wPost(:,j));
        
        for l = 1 : S
            xEstM(:,j,k) = xEstM(:,j,k) + xPost(:,l,j) * wPost(l,j); 
        end
        
        % state estimation
        xEst(:,k) = xEst(:,k) + xEstM(:,j,k);
        
        % mode specific state estimation
        xEstM(:,j,k) = xEstM(:,j,k) / mPost(j,k);
    end 
end

%% analysis of results

    % time vector 
    time = linspace(dt,T,t);

    % final mode estimation
    mode = zeros(1,t);
    for k = 1 : t
        maxP = max(mPost(:,k));
        maxMode = find(mPost(:,k) == maxP);
        if size(maxP,1) > 1
            mode(k) = mode(k-1);
        elseif (length(maxMode)) > 1
            mode(k) = mode(k-1);
        else
          mode(k) = (find(mPost(:,k) == maxP));  
        end
    end
    
    %RMS error of state estimation
    xRMS = (x - xEst).^2;

%% save the results
% A) in mat file
    %save('results.mat','time','mPost','x','xEstM','xEst','xRMS','mode') 

% B) in txt file for pgf plots
    % resultTable = [time; mode; x; xEst; xRMS];
    % fileID = fopen('results.txt','w');
    % fprintf(fileID,'%3s %1s %10s %10s %10s %10s %10s %10s %12s %12s %12s \n ','t','m', 'x1','x2','x3','xEst1','xEst2','xEst3','xRMS1', 'xRMS2', 'xRMS3');
    % fprintf(fileID,'%3.1f %1f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %12.10f %12.10f %12.10f \n',resultTable);
    % fclose(fileID);
    
%% plot results
    
figure(1) %State estimation
for i = 1 : n
    plot(time,x(i,:),'DisplayName',['x_' num2str(i)])
    hold on;
    plot(time,xEst(i,:),'--','DisplayName',['estimate x_' num2str(i)])
end
legend;
xlabel('time [s]')
ylabel('water level [m]')
hold off;

figure(2) %mode probability
for j = 1 : q
    plot(time,mPost(j,:),'DisplayName',['P(m^{(' num2str(j) ')})'])
    hold on;
end
legend;
xlabel('time [s]')
ylabel('posterior Mode probability')
hold off;

figure(3) %mode estimation
plot (time,mode);
xlabel('time [s]')
ylabel('Mode')
ylim ([0.8 q+0.2])

figure(4) %RMS error
for i = 1 : n
    plot(time,xRMS(i,:),'DisplayName',['RMS x_' num2str(i)])
    hold on;
end
legend;
xlabel('time [s]')
ylabel('water level RMS error')
hold off;

