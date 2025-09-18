function [optimalfitness,optimalSolution,RecgBestFitness] = WSAO(Particles_no,Max_iteration,LB,UB,dim,fobj)
%% WSAO —— Init
% Parameter Setting
nParticles = Particles_no; % num of agents
nDim = dim; % dim of agents
maxIterations = Max_iteration; % max Iteration
nNest= 5;  %num of ant Nests

R_k=0.4; %Fill Rate of n-dim Nests
Rnest = power(R_k,1/nDim)*(UB(1)-LB(1))*power(gamma(nDim/2+1),1/nDim)/(power(nNest,1/nDim)*sqrt(pi));  % Radius of Nest --Eq.5 

wMax = 0.9;
wMin = 0.1;
FADs=0.2;
% Adjust the lower and upper bounds lb and ub based on the input 
% (if the input is only one number, it is assumed that the lower bound is the same for each dimension).
if length(UB)== 1
    ub = ones(1,dim)*UB;
    lb = ones(1,dim)*LB;
else
    ub=UB;  %up bound 
    lb=LB;  %low bound
end

Vmax = 0.1*(ub(1)-lb(1));               % max velocity
Vmin = -Vmax;                      % min velocity

Xmin=repmat(ones(1,nDim).*lb,nParticles,1);  % low bound for size of nDim nParticles
Xmax=repmat(ones(1,nDim).*ub,nParticles,1);  % up bound for size of nDim nParticles

% Init the data record array
RecgBestFitness = zeros(1,maxIterations);
RecpBestFitness = zeros(maxIterations,nParticles);
RecgBestPositions = zeros(maxIterations,nDim);
RecpBestPositions = zeros(maxIterations,nParticles,nDim);
RecWeather = zeros(maxIterations,1);
RecWeathernum = zeros(3,1); % Record the number of each type of weather
RecAllSortedPoints = [];

% Init agents
positions = (ub-lb).*rand(nParticles, nDim)+lb; % init pos Eq.3
velocities = zeros(nParticles, nDim); % Init v
pBestPositions = positions; % init pBestpos
pBestFitness = inf(nParticles, 1); % init pBest
gBestPosition = zeros(1, nDim); % init gBestpos
gBestFitness = inf; % init gBest
fitness = inf(nParticles, 1);
tempfitness = fitness;
% Definition of Weather Conditions
SUNNY = 1;
CLOUDY = 2;
RAINY = 3;

% Init Weather State
weatherState = SUNNY;

% Check Bound and
% Calculate the fitness value 
for i = 1:nParticles
    positions(i,positions(i,:)>ub)=ub(1);
    positions(i,positions(i,:)<lb)=lb(1);
    fitness(i) = fobj(positions(i,:));
end

% Update the optimal position and fitness
improvedIndices = fitness < pBestFitness;
pBestPositions(improvedIndices, :) = positions(improvedIndices, :);
pBestFitness(improvedIndices) = fitness(improvedIndices);
[minFitness, minIndex] = min(fitness);
if minFitness < gBestFitness
    gBestPosition = positions(minIndex, :);
    gBestFitness = minFitness;
end

%% WSAO MianLoop
for iter = 1:maxIterations
    % Record data and  locate nests.
    % Create and updata the best historical ranking.As a candidate for the
    % nest location ranking.
    currentGenPoints = [fitness,positions];
    [~, sortedIdx] = sort(currentGenPoints(:, 1)); 
    sortedPoints = currentGenPoints(sortedIdx, :); % sorted Current Generation's Points
    
    % Insert the current point into the historical best points using merge sort
    i = 1; % Pointer to sortedPoints
    j = 1; % Pointer to RecAllSortedPoints
    m = size(RecAllSortedPoints, 1);
    n = size(sortedPoints, 1);
    % Initialize the result array
    result = zeros(m + n, dim+1);

    % Merge sort insertion
    k = 1; % Pointer to the result array
    while i <= n && j <= m
        if sortedPoints(i, 1) < RecAllSortedPoints(j, 1)
            result(k, :) = sortedPoints(i, :);
            i = i + 1;

        else
            result(k, :) = RecAllSortedPoints(j, :);
            j = j + 1;
        end
        k = k + 1;
    end
    % Deal with the remaining points
    while i <= n
        result(k, :) = sortedPoints(i, :);
        i = i + 1;
        k = k + 1;
    end
    
    while j <= m
        result(k, :) = RecAllSortedPoints(j, :);
        j = j + 1;
        k = k + 1;
    end
    RecAllSortedPoints = result(1:k-1, :);
    
    % Locate a new nest.    --Algorithm 1
    nests_now = [];
    % Select num of Nnests points from RecAllSortedPoints.
    for i = 1:size(RecAllSortedPoints, 1)
        canAdd = true;
        for j = 1:size(nests_now, 1)
            if euclideanDistance(RecAllSortedPoints(i, :), nests_now(j, :)) < Rnest
                canAdd = false;
                break;
            end
        end
        if canAdd
            nests_now = [nests_now; RecAllSortedPoints(i, :)];
        end
        if size(nests_now, 1) == nNest
            break;
        end
        % Replenish the nest
        if i == size(RecAllSortedPoints, 1)
            lacknum = nNest-size(nests_now, 1);
            suppleLords = RecAllSortedPoints(randi([1,size(RecAllSortedPoints, 1)],lacknum,1),:);
            nests_now = [nests_now;suppleLords];
        end
    end
    if iter == 1
        Nest_allocation = nests_now(randi([1,nNest],nParticles,1),2:end); %Nest allocation--Eq.8
    end
    % Records
    RecgBestFitness(iter) = gBestFitness;
    RecpBestFitness(iter,:) = pBestFitness;
    RecgBestPositions(iter,:) = gBestPosition;
    RecpBestPositions(iter,:,:) = positions;
    RecWeather(iter) = weatherState;
    RecWeathernum(weatherState) = RecWeathernum(weatherState) + 1;

    % WSAO agents updata 
    switch weatherState
        case SUNNY       %sunny state----Eq.7  
            Ws = wMax - iter .* ((wMax - wMin) / maxIterations);
            Cs1 = 1.6;
            Cs2 = 0.4+(1.6-0.4)*iter/maxIterations;
            velocities = Ws * velocities + ...
                Cs1 * repmat(rand(nParticles, 1),1,nDim) .* (pBestPositions - positions) + ...
                Cs2 * repmat(rand(nParticles, 1),1,nDim) .* (repmat(gBestPosition, nParticles, 1) - positions);
            velocities(velocities>Vmax) = Vmax;
            velocities(velocities<Vmin) = Vmin;
            positions = positions + velocities;
            

        case CLOUDY  %cloudy state----Eq.8  
            Nest_allocation = nests_now(randi([1,nNest],nParticles,1),2:end); %Nest allocation
            
            Wc= 0.1;
            Cc1 = 0.6;
            Cc2 = 2.0;
            velocities = Wc * velocities + ...
                Cc1 * rand(nParticles, nDim) .* (pBestPositions - positions) + ...
                Cc2 * rand(nParticles, nDim) .* (Nest_allocation - positions);
            velocities(velocities>Vmax) = Vmax;
            velocities(velocities<Vmin) = Vmin;
            positions = positions + velocities;

        case RAINY  % rainy state----Eq.9  

            Cr1 =1.6;
            Wr = 0.4;
            velocities = Wr * velocities + ...
                    Cr1 * rand(nParticles, nDim) .* (pBestPositions - Nest_allocation);
            velocities(velocities>Vmax) = Vmax;
            velocities(velocities<Vmin) = Vmin;
            positions = Nest_allocation + 0.5 * velocities;    

    end

    % The weather changes according to the state machine 
    weatherState = changeWeatherState(weatherState);  % --Eq.2

    
    % Check Bound and
    % Calculate the fitness value 
    for i = 1:nParticles
        positions(i,positions(i,:)>ub)=ub(1);
        positions(i,positions(i,:)<lb)=lb(1);
        fitness(i) = fobj(positions(i,:));
    end
    improvedIndices = fitness < pBestFitness;
    pBestPositions(improvedIndices, :) = positions(improvedIndices, :);
    pBestFitness(improvedIndices) = fitness(improvedIndices);

    %levy flight --Eq.11
    if rand()<FADs
%         U=rand(nParticles,nDim)<FADs;          
        cl1=(1-iter/maxIterations)^(2*iter/maxIterations);
%         temppositions=pBestPositions+cl1*((Xmin+rand(nParticles,nDim).*(Xmax-Xmin)).*U);
        temppositions=pBestPositions+cl1*((Xmin+rand(nParticles,nDim).*(Xmax-Xmin)));
    else
        r=rand();  Rs=size(positions,1);    cl2=0.2;
        stepsize=(cl2*(1-r)+r)*(pBestPositions(randperm(Rs),:)-pBestPositions(randperm(Rs),:));
        temppositions=pBestPositions+stepsize;
    end

    % Check Bound and
    % Calculate the fitness value 
    for i = 1:nParticles
        temppositions(i,positions(i,:)>ub)=ub(1);
        temppositions(i,positions(i,:)<lb)=lb(1);
        tempfitness(i) = fobj(temppositions(i,:));
    end
    improvedIndices = tempfitness < pBestFitness;
    positions(improvedIndices, :) = temppositions(improvedIndices, :);
    pBestPositions(improvedIndices, :) = temppositions(improvedIndices, :);
    pBestFitness(improvedIndices) = tempfitness(improvedIndices);
    fitness(improvedIndices) = tempfitness(improvedIndices);
    
    % Update the global optimal position and fitness
    [minFitness, minIndex] = min(fitness);
    if minFitness < gBestFitness
        gBestPosition = positions(minIndex, :);
        gBestFitness = minFitness;
    end
end

% Return the optimal value.
optimalSolution = gBestPosition;
optimalfitness = gBestFitness;
end

%% auxiliary functions
% auxiliary function：Weather Change
function newState = changeWeatherState(currentState)
    SUNNY = 1;
    CLOUDY = 2;
    RAINY = 3;

    s2s=0.6;
    s2c=0.25;
    s2r=0.15;
    c2s=0.25;
    c2c=0.25;
    c2r=0.5;
    r2s=0.25;
    r2c=0.05;
    r2r=0.7;
    
    rp=rand();
    switch currentState
        case SUNNY
            if rp<s2s
                newState=SUNNY;
            elseif rp<s2s+s2c
                newState=CLOUDY;
            else
                newState=RAINY;
            end
        case CLOUDY
            if rp<c2s
                newState=SUNNY;
            elseif rp<c2s+c2c
                newState=CLOUDY;
            else
                newState=RAINY;
            end
        case RAINY
            if rp<r2s
                newState=SUNNY;
            elseif rp<r2s+r2c
                newState=CLOUDY;
            else
                newState=RAINY;
            end
    end
end

% auxiliary function：euclideanDistance
function d = euclideanDistance(p1, p2)
    d = sqrt(sum((p1-p2).^2));
end

