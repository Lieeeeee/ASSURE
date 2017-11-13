function [mask, pntsArr] = drawLines(imgSize, P1, P2)
%generates a mask with lines, determine by P1 and P2 points

mask = zeros(imgSize);

P1 = double(P1);
P2 = double(P2);
pntsArr = generateLinesPoints(P1,P2);

for ii=1:size(P1,1)    
    % Determine x and y locations along the line
    xvalues = pntsArr(:,2,ii);
    yvalues = pntsArr(:,1,ii);
    % Replace the relevant values within the mask
    mask(sub2ind(size(mask), yvalues, xvalues)) = 1;

end

end