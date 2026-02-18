function throwError(logID,msg)
    % Log the error message
    fprintf(logID, '%s\n', msg); 
    % Throw the error to halt execution
    error(msg);
end