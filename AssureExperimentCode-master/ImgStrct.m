classdef ImgStrct
    %IMGSTRCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        img,
        imgFileName,
        frame,
        masks,
        masksFilesNames,
        uncMask,
        uncMasksFileName,
        cropBoundaries,
        dimensions
    end
    
    methods(Static)
        
        function res = getFrameImgStrct(imgStrct, frame)
            res = imgStrct;
            res.img = res.img(:,:,frame); 
            for z=1:(length(res.masks))
               res.masks{z} = res.masks{z}(:,:,frame); 
            end         
        end
        
        function res = cropImgStrct(imgStrct, xLeft,xRight,yTop,yBottom)
            res = imgStrct;
            res.img = imgStrct.img(xLeft:xRight,yTop:yBottom,:); 
            for z=1:(length(res.masks))
               res.masks{z} = res.masks{z}(xLeft:xRight,yTop:yBottom,:); 
            end       
        end
        
        function res = cropToSquare(imgStrct)
            res = imgStrct;
            [h,w] = size(res.img);
            newDim =  min(h,w);
            w = newDim;
            h = newDim;
            
            res.img = res.img(1:w,1:h,:); 
            for z=1:(length(res.masks))
               res.masks{z} = res.masks{z}(1:w,1:h,:); 
            end         
        end
        
        function strctOut = removeNovices(imgStrcts)
            strctOut = imgStrcts;
            for z=1:length(strctOut)
                masksToRemove = [];
                
                for k=1:length(strctOut{z}.masksFilesNames)
                    name = SUREExperimentCls.getRadiologistName(strctOut{z}.masksFilesNames{k});
                    seniority = SUREExperimentCls.getRadiologistSeniority(name);
                    if seniority ==1
                        masksToRemove = [masksToRemove,k];
                    end
                end
                strctOut{z}.masksFilesNames(masksToRemove) = [];
                strctOut{z}.masks(masksToRemove) = [];
            end
        end

        function strctOut = removeRadiologists(imgStrcts,radNames)
            strctOut = imgStrcts;
            for z=1:length(strctOut)
                masksToRemove = [];
                
                for k=1:length(strctOut{z}.masksFilesNames)
                    name = SUREExperimentCls.getRadiologistName(strctOut{z}.masksFilesNames{k});
                    removeCurrent = false;

                    for kk=1:length(radNames)
                        if isequal(name,radNames{kk})
                            removeCurrent = true;
                            break;
                        end
                    end
                    if removeCurrent
                        masksToRemove = [masksToRemove,k];
                    end
                end
                strctOut{z}.masksFilesNames(masksToRemove) = [];
                strctOut{z}.masks(masksToRemove) = [];
            end
        end
        function printAnnotatorsSt(imgStrcts)
            for z=1:length(imgStrcts)
                imgStrcts{z}.printAnnotators;
            end
        end
        
        function [ imgStrctCell ] = dbFileStrctToImgStrctCell( dbFiles, slices , is3DMode)
            %DBFILESTRCTTOIMGSTRCT Summary of this function goes here
            %   Detailed explanation goes here
            imgStrctCell = {};
            for t=1:length(dbFiles)
                %[t,length(dbFiles)]
                if isempty(slices{t})
                    continue;
                end
                tempStrcts = ImgStrct.generateImgStrct(dbFiles{t}.im,dbFiles{t}.segs,true,slices{t}, dbFiles{t}.windowing, is3DMode);
                if length(tempStrcts)>1
                    imgStrctCell = {imgStrctCell{:},tempStrcts{:}};
                else
                    imgStrctCell = {imgStrctCell{:},tempStrcts};
                end
            end
            
        end
        
        
        function displayImgCell(imgStrctCell,folderName)
            if ~exist(folderName,'dir')
                mkdir(folderName);
            end
            ii = randi(100);
            for t=1:length(imgStrctCell)
                ImgStrct.displayToFile(imgStrctCell{t},[folderName '\' num2str(t) '_' num2str(ii) '.png'])
            end
        end
        
        function displayToFile(imgStrct,fileName)
            
            [~,~,uncertaintyMask] = Utils.calcUnionIntersection(imgStrct.masks);
            close all;
            
            %segmentations
            LineDisplay.displaySegsOverlay(imgStrct.img,imgStrct.masks);
            im0 = LineDisplay.getCroppedFrameFromFigure();
            close all
            
            %uncertaintyMask
            LineDisplay.displaySegsOverlay(imgStrct.img,[],uncertaintyMask);
            im1 = LineDisplay.getCroppedFrameFromFigure();
            close all
            
            %sosna
            if ~isempty(imgStrct.uncMask)
                LineDisplay.displaySegsOverlay(imgStrct.img,[],imgStrct.uncMask);
                im2 = LineDisplay.getCroppedFrameFromFigure();
                close all
            else
                im2 = uint8(zeros(size(im0,1),1,3));
            end
            
            LineDisplay.displaySegsOverlay(imgStrct.img,zeros(size(imgStrct.img)));
            resizedIm = LineDisplay.getCroppedFrameFromFigure();
            close all
            combIm = [resizedIm,im0,im1,im2];
            
            imwrite(combIm,fileName);
            
        end
        
        function imgStrctCell = setContrast(imgStrctCell,contrastMap)
            for z = 1:length(imgStrctCell)
                [~,name,ext] = fileparts(imgStrctCell{z}.imgFileName);
                try
                    ran = contrastMap([name,ext]);
                    imgStrctCell{z}.img = mat2gray(imgStrctCell{z}.img,ran);
                catch
                    error(['couldnt find range: ' [name,ext]]);
                end
            end
            
        end
        
        function imgStrctOut = keepNMasks(imgStrct,N)
            imgStrctOut = imgStrct;
            imgStrctOut.masks(N+1:end) = [];
        end
        
        function imgStrcts = generateImgStrct(imName,segNames,cropFlag, relevantSlices, windowing, is3DMode)
            
            I = IO.loadFile(imName);
            if exist('windowing','var')
                I.img = mat2gray(I.img,windowing);
            end
            if ~exist('is3DMode','var')
                is3DMode = false;
            end
            segs = cell(length(segNames),1);
            
            %extracts dimensions
            if ~isequal(size(I.img),I.hdr.dime.dim(2:4))
                error 'cant extract dimensions'
            end
            dimensionsVal = I.hdr.dime.pixdim(2:4);
            
            calculateRelevantSlices = false;
            if ~exist('relevantSlices','var')
                calculateRelevantSlices = true;
                relevantSlices = 1:size(I.img,3);
            end
            allSegs = false(size(I.img));
            for t=1:length(segs)
                segs{t} = IO.loadFile(segNames{t});
                allSegs = allSegs | segs{t}.img>0;
            end
            
            if ~calculateRelevantSlices
                %ignores segs without all these slices
                tempSegs = {};
                tempSegNames = {};
                for t=1:length(segs)
                    if isequal(intersect(find(sum(sum(segs{t}.img,1),2)),relevantSlices),relevantSlices) || ...
                            isequal(intersect(find(sum(sum(segs{t}.img,1),2)),relevantSlices),relevantSlices')
                        tempSegs = {tempSegs{:},segs{t}};
                        tempSegNames = {tempSegNames{:},segNames{t}};
                    else
                        segNames{t}
                    end
                end
                segs = tempSegs;
                segNames = tempSegNames;
            end
            
            if cropFlag
                [minY,maxY,minX,maxX,minZ,maxZ] = Utils.extractMinVals(allSegs,25);
                crop = IO.extractCropOptStrctFromMinMaxVals(size(allSegs),minX,minY,minZ,maxX,maxY,maxZ);
                %crop strct determines how much we need to cut from each
                %dimension. we transform it to coordintaes
                crop.L = crop.cut_from_L+1;
                crop.R = size(allSegs,1)-crop.cut_from_R;
                crop.P = crop.cut_from_P+1;
                crop.A = size(allSegs,2)-crop.cut_from_A;
                crop.I = crop.cut_from_I+1;
                crop.S = size(allSegs,3)-crop.cut_from_S;
                if ~calculateRelevantSlices
                    %if we use relevant slices, we dont need to crop
                    crop.S = size(allSegs,3);
                    crop.I = 1;
                end
                I.img = I.img(crop.L:crop.R,crop.P:crop.A,crop.I:crop.S);
                for t=1:length(segs)
                    segs{t}.img = segs{t}.img(crop.L:crop.R,crop.P:crop.A,crop.I:crop.S)>0;
                end
                
            end
            
            if calculateRelevantSlices
                for t=1:length(segs)
                    relevantSlices = intersect(relevantSlices,find(sum(sum(segs{t}.img,1),2)>0));
                end
            end
            
            
            
            
            imgStrcts = cell(length(relevantSlices),1);
            if ~is3DMode
                for t=1:length(relevantSlices)
                    imgStrcts{t} = ImgStrct;
                    imgStrcts{t}.cropBoundaries = crop;
                    imgStrcts{t}.img = I.img(:,:,relevantSlices(t));
                    imgStrcts{t}.imgFileName = imName;
                    imgStrcts{t}.frame = relevantSlices(t);
                    imgStrcts{t}.dimensions = dimensionsVal;
                    for l=1:length(segs)
                        imgStrcts{t} = imgStrcts{t}.addSeg(segs{l}.img(:,:,imgStrcts{t}.frame)>0,segNames{l});
                    end
                end
            else
                imgStrcts = ImgStrct;
                imgStrcts.cropBoundaries = crop;
                imgStrcts.img = I.img(:,:,relevantSlices);
                imgStrcts.imgFileName = imName;
                imgStrcts.frame = relevantSlices;
                imgStrcts.dimensions = dimensionsVal;
                for l=1:length(segs)
                    imgStrcts = imgStrcts.addSeg(segs{l}.img(:,:,imgStrcts.frame)>0,segNames{l});
                end
            end
            
            
            
        end
        
        
    end
    methods
        function displayImgStrct(imgStrct)
            LineDisplay.displaySegsOverlay(imgStrct.img,imgStrct.masks);
        end
        function imgStrct = addSeg(imgStrct,seg,segName)
            imgStrct.masks{length(imgStrct.masks)+1} = seg;
            imgStrct.masksFilesNames{length(imgStrct.masksFilesNames)+1} = segName;
        end
        
        function imgStrct = addUncMask(imgStrct,uncMask,uncMaskName)
            imgStrct.uncMask = uncMask;
            imgStrct.uncMasksFileName = uncMaskName;
        end
        
        function printAnnotators(imgStrct)
            for z=1:length(imgStrct.masksFilesNames)
                [~,f] = fileparts(imgStrct.masksFilesNames{z});
                fprintf('%s \n',f)
                
            end
            fprintf('\n')
        end
    end
    
end

