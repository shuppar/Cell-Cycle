%% Function for generating maximum projections of stacks of tifs

function [stack, maxproj] = MaxNStack(fileFolder,dirProt, m, n)

    if nargin == 2
    
        m = 1; n = length(dirProt);
    
    else
        
    end

fileProt = {dirProt.name}';

stack = [];

for i = m:n
    
    stack = cat(3, stack, imread(fullfile(fileFolder, fileProt{i})));
    
end

maxproj = max(stack,[], 3);