function [filename] = get_filenameext(originalfilename,load_filename,dstcode)
% get filename extension: dataset identifier

% if district
if ~isnan(dstcode)
    filename = originalfilename + " district " + string(dstcode);
else
    filename = originalfilename;
    
    % if irrigation or rainfed or central
    pat = ["irrigationdst","rainfeddst"];

    if contains(load_filename,pat(1))
        filename = filename + " " + pat(1);
    elseif contains(load_filename,pat(2))
        filename = filename + " " + pat(2);
    else
        filename = filename + " all-districts" ;
    end
end


end