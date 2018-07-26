function experimentParams = readExperimentSpreadsheet(strain,includeAll)
    if nargin < 1
        strain = 'Ntsr1';
    end

    experimentParams = readtable('all_experiments.xlsx','Sheet',sprintf('%s maps',strain)); % TODO : top dir
    
    textFormattedNumericFields = {'Midline' 'ExcludeColumns' 'ExcludeRows' 'ExcludePixels' 'ProbeOrder'};
    
    for ii = 1:numel(textFormattedNumericFields)
        field = textFormattedNumericFields{ii};
        experimentParams.(field) = cellfun(@str2num,experimentParams.(field),'UniformOutput',false);
    end
    
    if nargin > 1 && includeAll
        return
    end
    
    experimentParams = experimentParams(strcmp(experimentParams.Include,'yes'),:);
end