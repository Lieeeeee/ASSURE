classdef ShapeAtlas
    %SHAPEATLAS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fileStrct;
        ignorePercentage = 0.1;
        latestRegisteredAtlas;
        latestImNii;
        %boundingBoxSegs;
    end
    
    methods (Static)
          
        function atlas = LeftCloseDefaultAtlas()
            opened = [2,4,7,10,13,14,17,19,20];
            extremelyOpen = [12];
            closed = [1,5,6, 8, 9,11,15,16,18];
            %closed = [1,5,6];
            invalid = [3];
            gtStrctCell = importdata('matStructs\kidneyLeftGTMovedCropped.mat');
            %gtStrctCell = importdata('matStructs\kidneyLeftGT.mat');
            gtStrctCell = IO.translateFileStrctToDifferentPC(gtStrctCell,'C:\Users\drorcohe\Desktop\thesisWork\','C:\Users\drorcohe\Desktop\SURECodeBase\SUREModel\');
            atlas = ShapeAtlas(gtStrctCell(closed));
            atlas = atlas.initAtlas();
        end
    end
    
    methods
      
        
        function [slices, slicesMask]= getSlicesToEvaluate(atlas, seg)
           slices = find(squeeze(sum(sum(seg,1),2))>0);
           nLayersToRemoveEachSide = floor((length(slices)*atlas.ignorePercentage)/2);
           slices(1:nLayersToRemoveEachSide) = [];
           slices(end-nLayersToRemoveEachSide+1) = [];
           
           slicesMask = seg;
           slicesMask(:,:,nLayersToRemoveEachSide) = 0;
        end
        
        function atlas = ShapeAtlas(fileStrct)
            atlas.fileStrct = fileStrct;
            %atlas.boundingBoxSegs = cell(size(fileStrct));
        end
        
        function atlasOut = initAtlas(atlas)
            atlasOut = ShapeAtlas(atlas.fileStrct);
            atlasOut.fileStrct = atlas.fileStrct;
            atlasOut.ignorePercentage = atlas.ignorePercentage;
            atlasOut.latestRegisteredAtlas = atlas.latestRegisteredAtlas;
            atlasOut.latestImNii = atlas.latestImNii;
            %for jj=1:length(atlas.boundingBoxSegs)
            %   seg = IO.loadFile(atlas.fileStrct{jj}.seg);
            %   seg.img = Utils.getBoundingBox(seg.img);
            %   seg.untouch=0;
            %   fileName = ['boundingBoxSeg' num2str(jj) '.nii.gz'];
            %   save_nii_gzip(seg,fileName);
            %   atlas.boundingBoxSegs{jj} = fileName;
            
            %end
        end
        
        function atlas = clearAtlas(atlas)
            %for jj=1:length(atlas.boundingBoxSegs)
            %    delete(atlas.boundingBoxSegs{jj});
            %end
            %atlas.boundingBoxSegs = [];
        end
        
        
        function [registeredAtlas, success, updatedAtlas] = registerToNii(atlas, imNiiFile, segNiiFile)
            if ~isempty(atlas.latestRegisteredAtlas) && isequal(atlas.latestImNii,imNiiFile)
                updatedAtlas = atlas;
                registeredAtlas = atlas.latestRegisteredAtlas;
                success = true;
                return;
            end
            
            
            TRANSFORMATION_NAME = 'temp\cpp.nii';
            IM_OUT_NAME = 'temp\outImNii.nii.gz';
            SEG_OUT_NAME = 'temp\outSegNii.nii.gz';
            SEG_BBOX_NAME =  'temp\segNiiBoundingBox.nii.gz';
            SEG_BBOX_NAME_ATLAS = 'temp\segNiiBoundingBoxAtlas.nii.gz';
            registeredSegs = cell(size(atlas.fileStrct));
            
            %array for evaluating result
            %registeredJaccard = zeros(size(atlas.fileStrct));
            registeredMI = zeros(size(atlas.fileStrct));
            transformationSize = zeros(size(atlas.fileStrct));
            
            segOrig = IO.loadFile(segNiiFile);
            segOrigMask = segOrig.img > 0;
            niiOrig = IO.loadFile(imNiiFile);
            
            %generates bounding boxes
            segBBoxIn = segOrig;
            segBBoxIn.untouch=0;
            segBBoxIn.img = Utils.getBoundingBox(segBBoxIn.img);
            %for z=1:size(segBBoxIn.img,3)
            %    segBBoxIn.img(:,:,z) = imdilate(segBBoxIn.img(:,:,z),strel('disk',4)); 
            %end
            save_nii_gzip(segBBoxIn,SEG_BBOX_NAME);
            
            
            
            for ii=1:length(atlas.fileStrct)
                ii
                
                bboxSeg = IO.loadFile(atlas.fileStrct{ii}.seg);
                bboxSeg.untouch=0;
                bboxSeg.img = Utils.getBoundingBox(bboxSeg.img);
                save_nii_gzip(segBBoxIn,SEG_BBOX_NAME_ATLAS);
                
                NiftyReg.registerNiftysAndMasks(imNiiFile,atlas.fileStrct{ii}.file,SEG_BBOX_NAME, SEG_BBOX_NAME_ATLAS,IM_OUT_NAME,SEG_OUT_NAME,TRANSFORMATION_NAME);
                NiftyReg.registerNiftysFromTransformation(segNiiFile,atlas.fileStrct{ii}.seg,SEG_OUT_NAME,TRANSFORMATION_NAME);
                
                registeredNii = IO.loadFile(IM_OUT_NAME);
                registeredSeg = IO.loadFile(SEG_OUT_NAME);
                closedMask = imopen(registeredSeg.img>0,strel('disk',1));
                for z =1:size(closedMask,3)
                    closedMask(:,:,z) = Utils.keepBiggestCC(closedMask(:,:,z));
                    closedMask(:,:,z) = imfill(closedMask(:,:,z),'holes');
                end
                figure,IO.imshow3D(closedMask);
                registeredSegs{ii} = closedMask;
                
                
                [~,targetRoi] = atlas.getSlicesToEvaluate(registeredSegs{ii});
                registeredMI(ii) = Utils.mutualInformation(registeredNii.img(targetRoi),niiOrig.img(targetRoi));
                
                [~,srcRoi] = atlas.getSlicesToEvaluate(bboxSeg.img);
                T = load_untouch_nii(TRANSFORMATION_NAME);
                T = imresize3d(squeeze(T.img),size(bboxSeg.img));
                %srcRoiMask3D = cat(4,srcRoi,cat(4,srcRoi,srcRoi));
                T = reshape(T,size(T,1)*size(T,2)*size(T,3),3);
                srcRoiReshaped = reshape(srcRoi,size(srcRoi,1)*size(srcRoi,2)*size(srcRoi,3),1);
                if size(srcRoiReshaped,1)~=size(T,1)
                    error 'err sizes'
                end
                T = T(srcRoiReshaped>0,:);
                transformationSize(ii) = mean(sqrt(sum(T.^2, 2)));
                
                registeredMI
                transformationSize
                delete TRANSFORMATION_NAME;
                %registeredJaccard(ii) = Unitest.jaccard(registeredSegs{ii},segOrigMask);
            end
            
            delete SEG_OUT_NAME;
            delete IM_OUT_NAME;
            delete TRANSFORMATION_NAME;
            delete BBOX_FILE_NAME;
            delete SEG_BBOX_NAME_ATLAS;
            
            normalizedMI = registeredMI/max(registeredMI(:));
            normalizedTSize = 1 - transformationSize/max(transformationSize(:));
            
            [~,maxInd] = max(registeredMI+normalizedTSize);
          %  [~,maxInd] = max(registeredJaccard);
            registeredAtlas = registeredSegs{maxInd};
            success = true;
            
            updatedAtlas = ShapeAtlas(atlas);
            updatedAtlas.latestRegisteredAtlas = registeredAtlas;
            updatedAtlas.latestImNii = imNiiFile;
            
            
        end
    end
    
end

