
% -*- texinfo -*- 
% @deftypefn {} {@var{protPositions} =} SimulateMultipleProtBinding_v4 (@var{DNAlength}, @var{minDist}, @var{maxDist}, @var(theta))
%  This function simulates protamines binding sequentially to a linear DNA molecule of length DNAlength
%  in an electrostatic potential of strength theta relative to the thermal
%  noise. It repeats the trial until there are numProt protamines within a distance maxDist of
%  one another, but no less than minDist away. It then returns their two
%  positions


% Author: rmcmi <rmcmi@LAPTOP-0665V2LG>
% Created: 2020-06-16

function [protPositions, boundPositions] = SimulateMultipleProtBinding_v6 (DNAlength, minDistance, maxDistance, theta, numProt)
  
  % The separation between the DNA and protamine when bound in nm
  protSeparation = 1;
  
  % The "width" of a protamine molecule bound to DNA
  l = 10;
  
  %repeatedly simulate numProt protamines binding
  for i = 1:DNAlength
     % Define the starting potential using the DNA alone
     potentialDimensionless = @(x) -1.*log(((1.-x).*DNAlength +sqrt((1.-x).^2*DNAlength.^2+protSeparation.^2))/(-x.*DNAlength+sqrt(x.^2.*DNAlength.^2+protSeparation.^2)));
     fnMaxXLoc = fminbnd(potentialDimensionless, 0, 1);
  
     % vector to store the fractional positions of protamines
     protPositions = [];
  
     % choose the first protPosition using the starting potential
     probDist = @(x) exp(-theta*potentialDimensionless(x));
     posOne = sampleDist(probDist,probDist(fnMaxXLoc),1,[0,1]);
     protPositions = [protPositions posOne];
    
     %vector to store the UNSORTED positions of protamines in the order
     %they bind
     boundPositions = [posOne];
     
     %Debug code to show the probability distribution as it updates
%    xVals = [0:0.001:1];
%    yVals = [];
%    for i = 1:length(xVals)
%      yVals(i) = probDist(xVals(i));  
%    end
%    figure;
%    plot(xVals, yVals);
     
     for n=2:numProt
         %Update the potential according to the last protamine that bound  
         potentialDimensionless = @(x) potentialDimensionless(x) + log(((protPositions(n-1)-x)*DNAlength+(l/2)+sqrt(((protPositions(n-1)-x)*DNAlength+(l/2)).^2+protSeparation.^2))/((protPositions(n-1)-x)*DNAlength-(l/2)+sqrt(((protPositions(n-1)-x)*DNAlength-(l/2)).^2+protSeparation.^2)));
         fnMaxXLoc = fminbnd(potentialDimensionless, 0, 1);

         % Update the probability distribution
         probDist = @(x) exp(-theta*potentialDimensionless(x));
  
         %Now sample from the updated distribution
         pos = sampleDist(probDist,probDist(fnMaxXLoc),1,[0,1]);
         protPositions = [protPositions pos];
         boundPositions = [boundPositions pos];
         
         
              %Debug code to show the probability distribution as it updates
%          xVals = [0:0.001:1];
%         yVals = [];
%         for i = 1:length(xVals)
%              yVals(i) = probDist(xVals(i));  
%         end
%         figure;
%          plot(xVals, yVals);
%          yVals;
     end
         % check if all protamines are within maxDistance but no less than
         % min distance
        protPositions = sort(protPositions);

        totDistance = (protPositions(numProt) - protPositions(1)) * DNAlength;
        
        %check that maxDistance condition is met
        if  totDistance < maxDistance
            1;
            
            % Now check that each pair of protamines is no less than
            % minDistance apart
            
            meetMinCondition = 1;
            for n=1:(length(boundPositions) - 1)
                if((protPositions(n+1) - protPositions(n)) * DNAlength) > minDistance
                    meetMinCondition = meetMinCondition + 1;
                end
            end
            if meetMinCondition == numProt
                return;
            end
        end 

 
    
  end 
end
