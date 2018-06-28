function plotPSTH(psth,isParameteric)
if isParametric
        [~,~,figureIndex] = unique([params.X params.Y],'rows');
        rowField = 'PulseWidth';
        colField = 'PercentPower';
    else
        [~,~,figureIndex] = unique([params.PulseWidth params.],'rows');
        rowField = 'Y';
        colField = 'X';
    end
    
    rowValues = unique(params.(rowField));
    colValues = unique(params.(colField));
end