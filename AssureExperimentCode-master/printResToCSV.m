function [  ] = printResToCSV( res, fileName )
%PRINTRESTOCSV print the result struct to a csv file using dlmwrite. every
%result is appended to this file as a single line
%   INPUT:
%       res - the results struct
%       fileName - the csv file to append the results to

    resToWrite = getResToPrint( res );
    dlmwrite(fileName, resToWrite, '-append', 'delimeter', '\t');
end

