function matlab_syntax = processJsonSyntax(filepath, opts)
arguments
    filepath        {isfile}        = './matlab-syntax.json'
    opts.saveFile   {islogical}     = true
end

% Decode json file
%   The file is typically a oniguruma YAML syntax file converted into a
%   json file.
mat_syntax_str = readlines(filepath);
matlab_syntax = jsondecode(strjoin(mat_syntax_str));

% All contexts are structs of action-value pairs with substep
% instructions
defStruct = struct("action",[],"value",[],"substep",struct([]));

% Process each entry in out.contexts
contexts = matlab_syntax.contexts;
contextNames = fieldnames(contexts);
variables = matlab_syntax.variables;

for idx = 1:length(contextNames)

    % Context name and values
    contextName = contextNames{idx};
    contextValues = contexts.(contextName);
    
    % Parse context values
    try 
        contextStruct = parseContext(contextValues);
    catch ME
        fprintf('Values for context %s is invalid', contextName);
        rethrow(ME)
    end

    % Substitute variables
    try
        contextStruct = subVariables(contextStruct);
    catch ME
        fprintf('Variable substitution for context %s valied', contextName);
        rethrow(ME)
    end
    
    % Store final contextStruct
    contexts.(contextName) = contextStruct;

end

% replace out.contexts
matlab_syntax.contexts = contexts;

% Save out in ../private
if opts.saveFile
    save(fullfile('..','private',"matlab-syntax.mat"),'matlab_syntax');
end


    function contextStruct = parseContext(contextValue)

        % Make a copy of the default struct
        contextStruct = defStruct;

        
        if isstruct(contextValue) && isscalar(contextValue)
            % scalar struct

            % Extract all actions
            actionNames = fieldnames(contextValue);
            for i = 1:length(actionNames)

                actionName = actionNames{i};

                switch actionName
                    case {'match', 'include', 'meta_scope', 'clear_scopes', ...
                            'meta_include_prototype', 'meta_content_scope'}
                        contextStruct.action = actionName;
                        contextStruct.value = contextValue.(actionName);
                        
                        
                        if actionName ~= "match" && ischar(contextStruct.value)
                            % Replace `-` with `_` for non match or char type
                            % values
                            contextStruct.value = strrep(contextStruct.value, '-', '_');
                        
                        elseif actionName == "match"
                            % Replace \b (word boundary), assuming they only
                            % appear at the ends of a expression
                            expr = contextStruct.value;
                            if startsWith(expr, '\b'), expr(1:2) = '\<'; end
                            if endsWith(expr, '\b'), expr(end-1:end) = '\>'; end
                            contextStruct.value = expr;

                        end
                    
                    otherwise
                        substepStruct = defStruct;
                        substepStruct.action = actionName;
                        switch actionName
                            case 'captures'
                                substepStruct.value = string(struct2cell(contextValue.captures));
                            otherwise
                                if ischar(contextValue.(actionName))
                                    substepStruct.value = strrep(contextValue.(actionName), '-', '_');
                                else
                                    substepStruct.value = contextValue.(actionName);
                                end
                        end
                        
                        % Append substepStruct
                        if isempty(contextStruct.substep)
                            contextStruct.substep =  substepStruct;
                        else
                            contextStruct.substep(end+1) =  substepStruct;
                        end
                end
            end


        elseif isstruct(contextValue)
            % non-scalar struct

            % Loop through each entry
            for i = length(contextValue):-1:1
                contextStruct(i) = parseContext(contextValue(i));
            end

        elseif iscell(contextValue)
            % cell array

            % Loop through each entry
            for i = length(contextValue):-1:1
                contextStruct(i) = parseContext(contextValue{i});
            end

        else
            error('Mismatch context value type')    
        end

    end

    function contextStruct = subVariables(contextStruct)

        % regex expression
        expr = '{{2}(?<var>[\w-]*)}{2}';

        % Iterate through contextStruct actions
        numActions = length(contextStruct);
        for actionIdx = 1:numActions
            
            % Retrieve name and value of action entry
            actionName = contextStruct(actionIdx).action;
            actionValue = contextStruct(actionIdx).value;

            if ischar(actionValue)
                
                % Match with regex
                [vars, valueParts] = regexp(actionValue, expr, "names","split","noemptymatch");
                % New actionValue
                actionValue_new = '';

                % Process each match
                for varIdx = 1:length(vars)
                    
                    % Retrieve variable value
                    varValue = variables.(vars(varIdx).var);

                    % Rebuild value using valueParts
                    actionValue_new = strcat(actionValue_new,valueParts{varIdx},varValue);
                end

                % Append last element of valueParts
                actionValue_new = strcat(actionValue_new,valueParts{end});

                % Store sub-ed value
                contextStruct(actionIdx).value = actionValue_new;

            elseif iscellstr(actionValue)
                %TODO: implement variable sub for cell arrays.

            end

        end

    end



end