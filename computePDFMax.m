function pdfMax = computePDFMax(currentDist,protPos)
%   Vector for storing the local minima
    localMinima = zeros(1,length(protPos)+1);
    
%   All basins of attraction for the local minima
    searchBounds = [0,sort(protPos),1];
    
%   Find the local minimum corresponding to each basin
    for i = 1:length(localMinima)
        [~,localMinima(i)] = fminbnd(@(x) -currentDist(x),searchBounds(i),searchBounds(i+1));
    end
    
%   Determine the pdf maximum
    pdfMax = max(-localMinima);
end