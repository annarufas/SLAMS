function [quot,err] = calculateUncertaintyInQuotient(xVals,yVals,exVals,eyVals)

    % This function computes the propagated uncertainty of a quotient for 1D 
    % or 2D input arrays of data values and their corresponding uncertainties, 
    % based on the general error propagation formula:
    % err = sqrt(ex^2*diff(q,x)^2 + ey^2*diff(q,y)^2)
    %
    % This function allows the input of 1D or 2D arrays, where rows are
    % the different values to be averaged for each column independently.

    % Check if all input arrays have the same size
    [nRows,nCols] = size(xVals);
    if any([size(yVals,1), size(yVals,2)] ~= [nRows, nCols]) || ...
       any([size(exVals,1), size(exVals,2)] ~= [nRows, nCols]) || ...
       any([size(eyVals,1), size(eyVals,2)] ~= [nRows, nCols])
        error('Input arrays must have the same dimensions');
    end

    % Initialise output arrays
    quot = NaN(nRows,nCols); 
    err = NaN(nRows,nCols); 
    
    % Identify where both xVals and yVals are not NaN and yVals are not
    % zero (makes the quotient infinite)
    validIndices = ~isnan(xVals) & ~isnan(yVals) & (yVals ~= 0);
    
    % Calculate the quotient for valid entries
    quot(validIndices) = xVals(validIndices)./yVals(validIndices);
    
    % Apply the error propagation formula for valid entries
    err(validIndices) = quot(validIndices) .*...
        sqrt((exVals(validIndices)./xVals(validIndices)).^2 + (eyVals(validIndices)./yVals(validIndices)).^2);

% Check out how the general error propagation formula works with the Matlab symbolic tool

% mA = 7.1;
% mB = 2.5;
% mC = 3.3;
% 
% eA = 3.9;
% eB = 2.0;
% eC = 1;
% 
% syms mA mB mC eA eB eC
% m = (mA + mB + mC)/3;
% e = sqrt(eA^2*diff(m,mA)^2+eB^2*diff(m,mB)^2+eC^2*diff(m,mC)^2); % general error formula
% avg = double(subs(m,{mA,mB,mC},{7.1,2.5,3.3}))
% err = double(subs(e,{mA,mB,mC,eA,eB,eC},{7.1,2.5,3.3,3.9,2.0,1.0})) % 1.6,2.9,3,0.88 

end