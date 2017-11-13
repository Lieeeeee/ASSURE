function [ pOut ] = simplifyPoly( p,nPoints )
%SIMPLIFYPOLY Summary of this function goes here
%   Detailed explanation goes here

T = 0;
jFactor = 1;
pOut = dpsimplify(p,T);
l = length(pOut);
if l < nPoints
    error('nPoints too big');
end
maxIters = 30;
c = 0;
while l~=nPoints && c<maxIters
    if l < nPoints
        jFactor = jFactor/2;
        T = T - jFactor;
    else
        T = T+jFactor;
    end
    pOut = dpsimplify(p,T);
    l = length(pOut);
    c = c + 1;
end

if c==maxIters
    if l > nPoints
        [pOut] = simplifyPoly(p,nPoints-1);
    end
    l = length(pOut);
    for ii=1:nPoints-l
        pOut = addPoint(pOut);
    end
end

if(length(pOut)~=nPoints)
    a = 3;
end

end

function pOut = addPoint(p)
    pOut = zeros(size(p,1)+1,size(p,2));
    pOut(1:end-1,:) = p;
	pOut(end,:) = (p(end,:)+p(1,:))*0.5; 
end