function e = calculateUncertaintyInAverage(xVals,eVals)

    % This function computes the propagated uncertainty of the mean for 1D 
    % or 2D input arrays of data values and their corresponding uncertainties, 
    % based on the general error propagation formula:
    % e = sqrt(exi^2*diff(m,xi)^2 + exj^2*diff(m,xj)^2 + ... + exn^2*diff(m,xn)^2)
    %
    % This function allows the input of 1D or 2D arrays, where rows are
    % the different values to be averaged for each column independently.

    % Normalise input to ensure column-wise processing
    if isrow(xVals)
        xVals = xVals';  % convert row vector to column vector
        eVals = eVals';
    end

    % Ensure inputs are the same size
    assert(isequal(size(xVals), size(eVals)), 'xVals and eVals must have the same dimensions');

    [nRows,nCols] = size(xVals);
    
    % Initialise the variance matrix
    variance = NaN(nRows,nCols);

    % Compute variance only where xVals are valid
    validMask = ~isnan(xVals);
    variance(validMask) = eVals(validMask).^2;

    % Corrected sample size for each column (ignoring NaNs)
    corrN = sum(validMask, 1);  % 1 value per column

    % Avoid division by zero
    corrN(corrN == 0) = NaN;

    % Propagate the uncertainty for each column
    e = sqrt(sum(variance, 1, 'omitnan')) ./ corrN;

    % Return as column if input was a column, otherwise as scalar for row vectors
    if nCols == 1
        e = e(:);
    end

% Below is a demonstration of how the error formula above and the general
% error propagation formula (written using the Matlab symbolic tool)
% produce the same results
%
% x = 6.942;
% y = 6.959;
% z = 6.88;
% 
% ex = 0.020;
% ey = 0.019;
% ez = 0.025;
% 
% N = 3;
% 
% err1 = sqrt(ex^2+ey^2+ez^2)/N; % 0.0124
% mean = 6.9270;
% 
% syms x y z ex ey ez
% m = (x + y + z)/N;
% generr = sqrt(ex^2*diff(m,x)^2+ey^2*diff(m,y)^2+ez^2*diff(m,z)^2); % general error formula
% err2 = double(subs(generr,{x,y,z,ex,ey,ez},{6.942,6.959,6.88,0.020,0.019,0.025}));
%
% This is useful to read:
% https://stats.stackexchange.com/questions/71419/average-over-two-variables-why-do-standard-error-of-mean-and-error-propagation

end