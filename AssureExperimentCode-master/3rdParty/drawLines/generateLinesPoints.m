function [pntsArr] = generateLinesPoints(P1, P2)
%generates a mask with lines, determine by P1 and P2 points
%pntsArr - 

P1 = double(P1);
P2 = double(P2);

for ii=1:size(P1,1)
    x1 = P1(ii,2); y1 = P1(ii,1);
    x2 = P2(ii,2); y2 = P2(ii,1);
    
    % Distance (in pixels) between the two endpoints
    nPoints = getNPointsBetweenP1P2(x1,y1,x2,y2);
    
    if ~exist('pntsArr','var')
       pntsArr = zeros(nPoints,2,size(P1,1));
    end
    
    if nPoints ~= size(pntsArr,1)
        error 'lines should have same length' 
    end
    % Determine x and y locations along the line
    xvalues = round(linspace(x1, x2, nPoints));
    yvalues = round(linspace(y1, y2, nPoints));
    pntsArr(:,:,ii) = [yvalues',xvalues'];
end

end