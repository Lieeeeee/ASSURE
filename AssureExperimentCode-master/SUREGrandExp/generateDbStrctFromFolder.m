function [ dbFiles ] = generateDbStrctFromFolder(annotationsFolders, volFolder, filePrefixes, windowings)
%GENERATEDBSTRCT Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(annotationsFolders)
    annotationsFolders = {annotationsFolders};
end
dbFiles = {};
for z=1:length(filePrefixes)
    segNames = {};
    for t=1:length(annotationsFolders)
        curNames = regexpdir(annotationsFolders{t}, [filePrefixes{z} ,'(.|_)*'], true);
        segNames = {segNames{:},curNames{:}};
    end
    volName = regexpdir(volFolder, [filePrefixes{z} ,'(.|_)*'], true);
    dbFiles{z}.im = volName{1};
    dbFiles{z}.segs = segNames;
    dbFiles{z}.windowing = [windowings{z}(1),windowings{z}(2)];
    
    
    
end




end
