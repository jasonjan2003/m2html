function s = mfileparse(mfile, mdirs, names, opts)
%MFILEPARSE Parsing of an M-file to obtain synopsis, help and references
%  S = MFILEPARSE(MFILE, MDIRS, NAMES, OPTIONS) parses the M-file MFILE looking
%  for synopsis (function), H1 line, subroutines and todo tags (if any).
%  It also fills in a boolean array indicating whether MFILE calls M-files 
%  defined by MDIRS (M-files directories) AND NAMES (M-file names).
%  The input OPTIONS comes from M2HTML: fields used are 'verbose', 'global'
%  and 'todo'.
%  Output S is a structure whose fields are:
%     o synopsis: char array (empty if MFILE is a script)
%     o h1line: short one-line description into the first help line
%     o subroutine: cell array of char containing subroutines synopsis
%     o hrefs: boolean array with hrefs(i) = 1 if MFILE calls mdirs{i}/names{i}
%     o todo: structure containing information about potential todo tags
%
%  See also M2HTML

% Copyright (C) 2003 Guillaume Flandin
arguments
    mfile
    mdirs
    names
    opts    {isstruct}      = struct('verbose', true, ...
                                     'globalHypertextLinks', false, ...
                                     'todo', false);
end


%- Delimiters used in strtok: some of them may be useless (% " .), removed '.'
strtok_delim = sprintf(' \t\n\r(){}[]<>+-*~!|\\@&/,:;="''%%');

%- Global Hypertext Links option
%  If false, hypertext links are done only among functions in the same
%  directory.
if opts.globalHypertextLinks
    hrefnames = names;
else
    indhref = find(strcmp(fileparts(mfile),mdirs));
    hrefnames = names(indhref);
end

%- Open for reading the M-file
fid = fopen(mfile,'r');
it = 0; % line number

%- Initialize Output
s = struct('synopsis',   '', ...
           'h1line',     '', ...
           'subroutine', {{}}, ...
           'hrefs',      sparse(1,length(names)), ...
           'todo',       struct('line',[],'comment',{{}}), ...
           'ismex',      zeros(size(mexexts)));

%- Initialize flag for synopsis cont ('...')
flagsynopcont = false;

%- Look for synopsis and H1 line
%  Help is the first set of contiguous comment lines in an m-file
%  The H1 line is a short one-line description into the first help line
while 1

    % Read new line
    tline = fgetl(fid);
    
    % Break if not a char (ex. EOF)
    if ~ischar(tline), break, end

    % Increment current line number
    it = it + 1;

    % Trim whitespaces
    tline = strtrim(tline);
    
    %- Synopsis line
    if startsWith(tline,'function')
        s.synopsis = tline;
        % Check if synopsis line continues (ends with `...`)
        if endsWith(tline,'...')
            flagsynopcont = true;
            s.synopsis = strtrim(s.synopsis(1:end-3));
        end

    %- H1 Line
    elseif startsWith(tline,'%')
        % allow for the help lines to be before the synopsis
        if isempty(s.h1line)
            % Remove '%' and white space in front of line
            s.h1line = strip(strip(tline,'left','%'),'left');
        end
        % Save only the comment line closest to the synopsis
        if ~isempty(s.synopsis), break, end

    %- Go through empty lines
    elseif isempty(tline)
        % Do nothing

    %- Found code, or continued lines of synopsis
    else
        % Append continued synopsis line
        if flagsynopcont
            % There are more continued lines
            if endsWith(tline, '...')
                s.synopsis = [s.synopsis strtrim(tline(1:end-3))];
            
            % No more continued lines
            else
                s.synopsis = [s.synopsis tline];
                flagsynopcont = false;
            end
        % Code found. Move on to next section.
        else
            break;
        end
    end
end

%- Compute cross-references and extract subroutines
%  hrefs(i) is 1 if mfile calls mfiles{i} and 0 otherwise
while ischar(tline)
    % Remove blanks at both ends
    tline = strtrim(tline);
    
    % Split code into meaningful chunks
    splitc = splitcode(tline);
    for j=1:length(splitc)

        if isempty(splitc{j}) || ...
            splitc{j}(1) == "'" || ...
            contains(splitc{j},'...')
            % Forget about empty lines, char strings or conts
            % Do nothing
        
        elseif splitc{j}(1) == '%'
            % Cross-references are not taken into account in comments
            % Just look for potential `% TODO` or `% FIXME` line
            if opts.todo
                todo_ptrn = "%\s*(TODO|FIXME)(?:\s|:)*(?<TODO>.*)";
                regexp_results = regexp(splitc{j}, todo_ptrn, "names");
                if ~isempty(regexp_results)
                    s.todo.line(end+1) = it;
                    s.todo.comment{end+1} = regexp_results.TODO;
                end
            end
        else
            % detect if this line is a declaration of a subroutine
            if startsWith(splitc{j}, 'function')
                s.subroutine{end+1} = splitc{j};
            else
                % get list of variables and functions
                symbol = {};
                while 1
                    [t,splitc{j}] = strtok(splitc{j},strtok_delim);
                    if isempty(t), break, end
                    symbol{end+1} = t;
                end
                if opts.globalHypertextLinks
                    s.hrefs = s.hrefs + ismember(hrefnames,symbol);
                else
                    if ~isempty(indhref)
                        s.hrefs(indhref) = s.hrefs(1,indhref) + ...
                                           ismember(hrefnames,symbol);
                    end
                end
            end
        end
    end
    tline = fgetl(fid);
    it = it + 1;
end	

fclose(fid);

%- Look for MEX files
[pathstr,name] = fileparts(mfile);
samename = dir(fullfile(pathstr,[name	'.*']));
samename = {samename.name};
ext = cell(1,length(samename));
for i=1:length(samename)
    [~, ~, ext{i}] = fileparts(samename{i});
    if numel(ext{i}) && ext{i}(1) == '.'
        ext{i}(1) = [];
    end
end
allexts = mexexts;
s.ismex = ismember({allexts.ext},ext);
