function [psth,sdf,mpsth,msdf,params] = loadPSTH(folder,isParametric)
    if nargin < 1
        if nargin == 1
            isParametric = folder;
        else
            isParametric = false;
        end
            
        folder = pwd;
    end
end