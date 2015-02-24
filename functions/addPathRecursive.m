% adds folder specified by path as well as all subfolders

function [] = addPathRecursive(path)

if exist(path, 'dir')
    
    % first add the original folder
    addpath(path);
    
    % get a list of everything in the folder
    dirList = dir(path);
    
    if length(dirList) > 2
    
        % first two will always be . and ..
        for i = 3:length(dirList)
            % if we have found a subfolder, add it
            if dirList(i,1).isdir
                addPathRecursive([path '/' dirList(i,1).name])
            end
        end
        
    end
    
else
    error('notDirectory', 'Input path must be a directory');
end