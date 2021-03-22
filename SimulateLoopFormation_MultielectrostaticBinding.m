function [ startSites ] = SimulateLoopFormation_MultielectrostaticBinding(N, a, b, d, t, theta, minDistance, maxDistance, numProt, numLoop )
% SIMULATELOOPFORMATION_MULTIELECTROSTATICBINDING Given a number of segments and 
% the parameters for a gamma distribution, this function computes an expected 
% distribution of loop start sites.
%   
%   This program conducts t trials in which it randomly
%   selects a loop start site from the N segments of the chain, but prefers to be in the middle.
%   The binding is set via an electrostatic model in which the DNA is a line charge
%   and the protamine is a point charge. We use Boltzmann statistics to convert
%   the potential into a probability. It simulates protamines binding to
%   DNA and only terminates when numProt of them bind within maxDistance, but no
%   less than minDistance. The protamine closest to any end is the center of
%   a loop, and we branch out half of the circumference in each direction.
    

%   Theta is the ratio of the magnitude of the electrostatic
%   energy to the thermal energy. We've been using something on the order of 1
%   You should set N to the length of the DNA in nm. a and b are alpha and
%   beta for the gamma distribution, while d is an x-offset.
%   maxDistance is the order of one persistence length. We're not currently using
%   minDistance, but left it in in case people want to experiment with it.
%   NumProt is the number of protamines to bind, and numLoop is the number
%   of loops. 

%   This code simulates the electrostatic binding model if numProt=1 and
%   the electrostatic multibinding model if numProt > 1. Note that if
%   numProt=1, then minDistance and maxDistance do not effect the results
%   at all.

startSites = zeros(1, t);
      
%Change this if the image is not 2 um x 2 um and 512 pixels x 512 pixels or 1 um x 1 um
%and 256 pixels x 256 pixels. Variable has units nm/pixel
pixelSize = 3.9;

%3 pixels times nm/pixel divided by nm gives a dimensionless binsize.
binSize = 3*pixelSize/N;

% Set bins for histogram
edges = [0:binSize:0.5];

%initialize rng
rng('default')

%conduct t trials
for i=1:t
    %randomly choose numProts positions.
    [~, bindingPositions] = (SimulateMultipleProtBinding_v6 (N, minDistance, maxDistance, theta, numProt));
    
    %randomly choose a number of steps
    numSteps = 0;
    for n=1:numLoop
        numSteps = numSteps + 3.14159* (gamrnd(a, 1/b) + d);
    end
  
    %determine which protamine is closest to the end
    bindingPositions = sort(bindingPositions);
    if bindingPositions(1) < (1-bindingPositions(numProt))
        midpoint = bindingPositions(1) * N;
    else
        midpoint = bindingPositions(numProt) * N;
    end
      
    leftEnd = midpoint - (numSteps / 2);
    rightEnd = midpoint + (numSteps / 2);
    
    %if endPos is off of the molecule, then we say that the start site is
    %0.
    distToRightEnd = N - rightEnd;
    if(leftEnd < 1 | rightEnd > N)
        startSites(i) = 0;
        continue;
        
    elseif(leftEnd < distToRightEnd)
        startSites(i) = leftEnd / N;
        
    elseif(leftEnd > distToRightEnd)
        startSites(i) = distToRightEnd / N;
        
    end

    
end 
startSites = startSites';
%make the histogram
figure('Name', ['DNA Loop Start Site Simulation with ' num2str(N) ' segments' ])
startSitesHist = histogram(startSites, edges); 
xlabel('Normalized Loop Start Site')
ylabel('N')
title(['DNA Loop Start Site Simulation with ' num2str(N) ' segments' ])
end

