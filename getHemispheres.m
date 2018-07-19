function [left,right] = getHemispheres(columns,midline)
    isNotBetweenPixels = all(round(midline) == midline);
    
    rightestLeft = floor(midline(1))-isNotBetweenPixels;
    leftestRight = ceil(midline(end))+isNotBetweenPixels;
    
    left = 1:rightestLeft;
    right = leftestRight:columns;
    
    nL = numel(left);
    nR = numel(right);
    
    if nL < nR
        right(end-(nR-nL)+1:end) = [];
    elseif nL > nR
        left(1:(nL-nR)) = [];
    end
end