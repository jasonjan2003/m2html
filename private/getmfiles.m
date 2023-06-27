function mfiles = getmfiles(mdirs, mfiles, recursive,ignoredDir)
% Extract M-files from a list of directories and/or M-files

if nargin < 4, ignoredDir = {}; end
for i=1:length(mdirs)
    currentdir = fullfile(pwd, mdirs{i});
    if exist(currentdir) == 2 % M-file
        mfiles{end+1} = mdirs{i};
    elseif exist(currentdir) == 7 % Directory
        d = dir(fullfile(currentdir, '*.m'));
        d = {d(~[d.isdir]).name};
        for j=1:length(d)
            %- don't take care of files containing ','
            %  probably a sccs file...
            if isempty(findstr(',',d{j}))
                mfiles{end+1} = fullfile(mdirs{i}, d{j});
            end
        end
        if recursive
            d = dir(currentdir);
            d = {d([d.isdir]).name};
            d = {d{~ismember(d,{'.' '..' ignoredDir{:}})}};
            for j=1:length(d)
                mfiles = getmfiles(cellstr(fullfile(mdirs{i},d{j})), ...
                                   mfiles, recursive);
            end
        end
    else
        fprintf('Warning: Unprocessed file %s.\n',mdirs{i});
        if ~isempty(strmatch('/',mdirs{i})) || findstr(':',mdirs{i})
            fprintf('         Use relative paths in ''mfiles'' option\n');
        end 
    end
end