function experimentParams = readExperimentSpreadsheet(includeAll)
    experimentParams = readtable('all_experiments.xlsx'); % TODO : top dir
    
    textFormattedNumericFields = {'Midline' 'ExcludeColumns' 'ExcludeRows' 'ExcludePixels'};
    
    for ii = 1:numel(textFormattedNumericFields)
        field = textFormattedNumericFields{ii};
        experimentParams.(field) = cellfun(@str2num,experimentParams.(field),'UniformOutput',false);
    end
    
    if nargin > 0 && includeAll
        return
    end
    
    experimentParams = experimentParams(strcmp(experimentParams.Include,'yes'),:);
end