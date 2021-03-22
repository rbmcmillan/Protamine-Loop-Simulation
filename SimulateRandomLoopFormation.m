function [ startSites ] = SimulateRandomLoopFormation(N, a, b, d, t, numLoop )
% SIMULATERANDOMLOOPFORMATION Given a number of segments and an average loop
% circumference, this function computes an expected distribution of loop
% start sites, assuming that loop formation is stochastic and unbiased.
%   
%   This program conducts t trials in which it randomly
%   selects a loop start site from the N segments of the chain. It then
%   chooses a loop size, which is drawn from a gamma distribution with
%   parameters a and b and x-offset d. We assume that the loop starts 180 degrees out of phase 
%   from where it ends, so we step a distance of half of the circumference to the right and to
%   the left. Finally, it computes the
%   distance between the nearest end of the DNA and the start site


%   This program simulates the random looping model.

startSites = [];
      
%Change this if the image is not 2 um x 2 um and 512 pixels x 512 pixels or 1 um x 1 um
%and 256 pixels x 256 pixels. Variable has units nm/pixel
pixelSize = 3.9;

%3 pixels times nm/pixel divided by nm gives a dimensionless binsize.
binSize = 3*pixelSize/N;

% Set bins for histogram
edges = [0:binSize:0.5];


%conduct t trials
for i=1:t
    %randomly choose a start position
    startPos = unidrnd(N);
    
    %randomly choose a number of steps
    numSteps = 0;
    for n=1:numLoop
        numSteps = numSteps + 3.14159* (gamrnd(a, 1/b) + d);
    end
    
    endPos = startPos + numSteps;
    
    leftEnd = startPos - (numSteps / 2);
    rightEnd = startPos + (numSteps / 2);
    
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
%make the histogram
figure('Name', ['DNA Loop Start Site Simulation with ' num2str(N) ' segments' ])
startSitesHist = histogram(startSites, edges); 
xlabel('Normalized Loop Start Site')
ylabel('N')
title(['DNA Loop Start Site Simulation with ' num2str(N) ' segments' ])

startSites = startSites';
end

