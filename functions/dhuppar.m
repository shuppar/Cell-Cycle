%% dhuppar function

function out = dhuppar(nrows, x)

if nrows - x > 0
    out = x;
else
    out = nrows;
    
end