%% shuppar function. It makes sure that the cropping region does not go beyond the image size.

function out = shuppar(x)

if x > 0
    out = x;
else
    out = 1;
    
end
