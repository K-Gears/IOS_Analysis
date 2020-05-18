function y = getRGBcomponent(x)
    % getRGBcomponent: convert fitting weight parameters to CData in a RGB
    %       image. Arterial component in red, and venous component in blue.
    %
    % USAGES:
    %       y = getRGBcomponent(x);
    %
    % INPUTS:
    %       x: fitting weight for arterial component or venous component
    %          N by 1 vector
    %
    % OUTPUTS:
    %       y: RGB component for arterial component or venous component
    %          N by 1 vector
    % 
    % Last Modified: 10-14-2015
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    
    y = x; % duplicate the input parameter and work on the copy

    y(x>0) = 0; % set all positive parameters to zero
    yMin = max(y(y~=0)); % get the max of all negative elements
    y = y-yMin; % correct baseline using the max of negative elements
    
    
%     % FIND 1 PERCENTILE OF ALL THE NEGATIVE ELEMENTS: METHOD 1 (Bing's method)
%     yNum = sum(y<0); % the number of all negative elements (i.e., remove baseline values)
%     [yCount, yBin] = hist(y(y<0),100); % histogram with 100 equal bins 
%     % Get the 1 percentile of these negative elements
%     for i = 1:100
%         if sum(yCount(1:i)) >= yNum/100
%             yMax = yBin(i-1);
%             break;
%         end
%     end

    
    % FIND 1 PERCENTILE OF ALL THE NEGATIVE ELEMENTS: METHOD 2 (Qingguang's method)    
    yMax = prctile(y(y<0),1); % find top 1 percentile
    
    y = y/yMax; % normalize
    y(y>1) = 1; % set element greater than 1 to 1
    y(y<0) = 0; % set element less than 0 to 0, this may set the arterial pixels in frontal cortex to 0
end