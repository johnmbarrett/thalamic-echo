function J = jet3(m)
% derived from Matlab's JET colormap
% this one doesn't have a black level
% 
% modified (gs)

%JET    Variant of HSV.
%   JET(M), a variant of HSV(M), is the colormap used with the
%   NCSA fluid jet image.
%   JET, by itself, is the same length as the current colormap.
%   Use COLORMAP(JET).
%
%   See also HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   C. B. Moler, 5-10-91, 8-19-92.
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.6 $  $Date: 2001/04/15 11:58:59 $
%




if nargin < 1, m = size(get(gcf,'colormap'),1); end
if m>15
    n = max(round(m/3.1),1);
    x = (1:n)'/n;
    y = (n/2:n)'/n;
    e = ones(length(x),1);
    r = [0*y; 0*e; x; e; flipud(y)];
    g = [0*y; x; e; flipud(x); 0*y];
    b = [y; e; flipud(x); 0*e; 0*y];
    J = [r g b];
    while size(J,1) > m
        J(1,:) = [];
        if size(J,1) > m, J(size(J,1),:) = []; end
    end
else
    if m==1
        J = [0 0 1
            1 0 0];
    end
    if m==2
        J = [0 0 1;
            .5 1 .5;
            1 0 0];
    end
    if m == 3
        J = [0 0 1;
            0 1 1;
            1 1 0;
            1 0 0];
    end            
end
