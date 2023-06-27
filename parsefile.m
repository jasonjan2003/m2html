function [outputArg1,outputArg2] = parsefile(code)
%PARSEFILE Summary of this function goes here
%   Detailed explanation goes here
arguments
    code {isstring} = readlines('../TwoPhaseSolver/+Inputs/Input.m')
end

defStruct = struct("action",[],"value",[],"substep",struct([]));
scopeStruct = struct("line",[],"start",[],"end",[],"scope",[],"meta_scope",[]);

% Load syntax file
load('private\matlab-syntax.mat','matlab_syntax');
contexts = matlab_syntax.contexts;

% Context stack with `main` as root
stack = "main";
scope = scopeStruct;
context = compileContext(stack(end));

% TODO: Figure out how to extract useful information

% Loop through each line of code
for lineIdx = 1:length(code)

    % Index line
    codeLine = code(lineIdx);

    % Parse line
    parseLine(codeLine);
    

end

    function compiledContext = compileContext(contextName)
    

        % Retrieve original context
        originalContext = contexts.(contextName);
        numActions = length(originalContext);

        % Setup compiledContext as empty array
        compiledContext = originalContext([]);

        % Flag for having `include` actions
        hasInclude = true;

        while 1

            % Loop through original context
            for actionIdx = 1:numActions
                
                if isempty(originalContext(actionIdx).action)
                    disp('pause')
                end
                switch originalContext(actionIdx).action
                    case 'include'
                        % "Include", only supports 1st-level includes
                        hasInclude = true;

                        % Retrieve context to include
                        subContext = contexts.(originalContext(actionIdx).value);
        
                        % length of subContext
                        numsubActions = length(subContext);
        
                        % Append subactions to compiledContext
                        compiledContext(end+(1:numsubActions)) = subContext;
    
                    otherwise
                        % simply append
                        compiledContext(end+1) = originalContext(actionIdx);

                end
    
            end

            if hasInclude
                % Reset compiledContext and hasInclude flag
                % Replace originalContext with compiledContext
                originalContext = compiledContext;
                compiledContext = originalContext([]);
                numActions = length(originalContext);
                hasInclude = false;
            else
                break
            end

            

        end

    end

    function parseLine(codeLine)
        
        % Search through line until line is exhausted
        while strlength(codeLine) > 1
         
            for actionIdx = 1:length(context)
                entry = context(actionIdx);
                meta_scope = [];
                meta_content_scope = [];

                switch entry.action
                    case 'match'
                        expr = entry.value;

                        % Attempty match
                        if expr == string('^(?:[A-Za-z]\w*)'), 
                            disp('pause')
                        end
                        [matchStartIdx,matchEndIdx, tokenExtents, codeLine] = regexp(codeLine,expr, "start", "end", "tokenExtents","split","lineanchors");
                        codeLine = strcat(codeLine{:});

                        % Continue if no matches are found
                        if isempty(matchStartIdx)
                            continue;
                        end
                        
                        % Process substeps
                        subSteps = entry.substep;
                        numSubSteps = length(subSteps);

                        for subStepIdx = 1:numSubSteps
                            
                            subStep = subSteps(subStepIdx);

                            switch subStep.action
                                case 'push'
                                    stack(end+1) = subStep.value;
                                case 'pop'
                                    for popIdx = 1:subStep.value
                                        %TODO: Here's a potential to pop
                                        %more than stack height
                                        stack(end) = [];
                                    end
                                case 'set'
                                    % Pop last on stack
                                    stack(end) = [];
                                    % Set context
                                    stack(end+1) = subStep.value;
                                case 'scope'
                                    newScope = scopeStruct;
                                    newScope.line = lineIdx;
                                    newScope.start = matchStartIdx;
                                    newScope.end = matchEndIdx;
                                    newScope.scope = subStep.value;
                                    newScope.meta_scope = meta_scope;
                                    scope(end+1) = newScope;
                                case 'captures'
                                    if length(matchStartIdx) ~= length(subStep.value)
                                        error('Mismatch captures array size: %u vs %u', ...
                                            length(matchStartIdx), ...
                                            length(subStep.value))
                                    end
                                    for capturesIdx = 1:length(subStep.value)
                                        
                                        % New scope for each match 
                                        newScope = scopeStruct;
                                        newScope.line = lineIdx;
                                        newScope.start = matchStartIdx(capturesIdx);
                                        newScope.end = matchEndIdx(capturesIdx);
                                        newScope.scope = subStep.value;
                                        newScope.meta_scope = meta_scope;
                                        scope(end+1) = newScope;
                                    end
                                otherwise
                                    error('New substep type: %s', subStep(subStepIdx).action)
                            end


                        end

                        % Done processing this line
                        context = compileContext(stack(end));
                        break;
                    case 'meta_scope'
                        meta_scope = entry.value;
                        continue;
                    case 'meta_content_scope'
                        %TODO: Not currently implemented
                        meta_content_scope = entry.value;
                        continue;
                    otherwise
                        error('New action type: %s', entry.action)
                end
            end
            % if actionIdx == length(context)
            % error('No matching action found: %s', stack(end))
        
        end
    end
    
        
end


        

