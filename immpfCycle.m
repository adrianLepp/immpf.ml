function [xPostKK,wPostKK] = immpfCycle(transP,xPostK,wPostK,y,systemSolver, systemParameter,q,S,n)
%immpf Cycle: IMMPF Algorithm for one time step k 
%homogen state space with n Dimensions, q Modes and S Particles
%Adrian Lepp, 25.02.2022

%___INPUTS:
%   transP              q x q Transition probability Matrix
%   xPostK              n x S x q posterior particle set for k-1
%   wPostK              S x q posterior weights for k-1
%   y                   measurement output for time step k
%   systemSovler        system Equations
%   systemParameter     parameters for system Equations: q dim struct
%   q                   number of Modes
%   S                   number of particles per Mode
%   n                   Dimension of state sptace vector

%___OUTPUTS:
%   xPostK              posterior particle set for k
%   wPostK              posterior weights for k
    
    % initialization
    mPrioKK = zeros(q,1);
    mPostK = mPrioKK;
    xPrioKK = zeros(n,S,q);
    xPostKK = xPrioKK;
    wPrioKK = zeros(S,q);
    wPostKK = wPrioKK;

    % posterior mode probability as sum of posterior weights
    for j = 1 : q
        mPostK(j) = sum(wPostK(:,j));
    end
    
    %% 1. mode switching
    for j = 1 : q 
        for i = 1 : q
            % prior mode probability
            mPrioKK(j) = mPrioKK(j) + mPostK(i)*transP(i,j); 
        end
    end
    
    % create particle set with S*q particles per Mode 
    % to approximate prior density after mode change
    xResampling = reshape(xPostK,[n,S*q]);
    wResampling = zeros(S*q,q);
    for j = 1 : q
        for i = 1 : q
            for l = 1 : S
                wResampling(l+(i-1)*S,j) = transP(i,j)*wPostK(l,i)  / mPrioKK(j) ;
            end
        end
    end
    % Reduce number of particles to S per Mode by Resampling
    for j = 1 : q
        if mPrioKK(j) <= q/(S*10) % if prior mode probabilty to low, take over posterior k-1
            xPrioKK(:,:,j) = xPostK(:,:,j);
        else
            for l = 1 : S
                pos = find(rand <= cumsum(wResampling(:,j)),1);
                xPrioKK(:,l,j) = xResampling(:,pos);       
            end
        end
        % new prior weights
        wPrioKK(:,j) = mPrioKK(j) ./ S;
    end 
       
    %% 2. Prediction
    for j = 1 : q
        for l = 1 : S
            % solve xKK = f(xK,uK,mKK) and yK = h(xKK,mKK) for alle particles
            [xPostKK(:,l,j),yTheor] = systemSolver(xPrioKK(:,l,j),systemParameter(j));
            
    %% 3. Correction       
            % calculate measurement probability p(y|x,m) by distance of
            % calculated yTheor and real measurement y with 
            % normal distribution
            pY = 1/((det(2*pi*systemParameter(j).sigmaY))^(0.5)) * exp(-0.5*(y - yTheor).' * inv(systemParameter(j).sigmaY)* (y - yTheor));
            
            % posterior weights
            wPostKK(l,j) = wPrioKK(l,j) * pY; 
            if wPostKK(l,j) < 1e-30
                wPostKK(l,j) = 1e-30;
            end
        end
    end
    % normalize weights
    summe = sum(wPostKK(:,:),'all');
    wPostKK(:,:) = wPostKK(:,:) ./summe;      
end