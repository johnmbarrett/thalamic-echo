function map = createMap(data,params)
    [rowValues,~,rowIndices] = unique(params.Y); % TODO : parameteric
    [colValues,~,colIndices] = unique(params.X);
    
    rows = numel(rowValues);
    cols = numel(colValues);
    
    assert(numel(data) == rows*cols,'Wrong number of data points in the map.');
    
    rowIndices = rows-rowIndices+1; % upside down
    
    map = zeros(rows,cols);
    map(sub2ind(size(map),rowIndices,colIndices)) = data;
end