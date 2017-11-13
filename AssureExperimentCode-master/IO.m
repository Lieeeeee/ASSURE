classdef IO
    %IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        LEFT_KIDNEY = 29663;
        RIGHT_KIDNEY = 29662;
        LIVER = 58;
        MAIN_GT_DIR = 'C:\Users\drorcohe\Dropbox\GroundTruth\';
        MAIN_ST_DIR = 'C:\Users\drorcohe\Dropbox\SilverCorpus\';
    end
    
    
    
    methods(Static)
        
        function saveAsTiff(tiffName, im)
            for K=1:size(im,3)
                currentSlice = squeeze(im(:,:,K,:));
                imwrite(currentSlice, tiffName, 'WriteMode', 'append');
            end
        end
        
        function uniteVideos(vid1,vid2,outName,FrameRate)
            if ~exist('FrameRate','var')
                FrameRate = 1;
            end
            v1 = VideoReader(vid1);
            v2 = VideoReader(vid2);
            v = VideoWriter(outName);
            v.FrameRate=FrameRate;
            v.open();
            while hasFrame(v1) && hasFrame(v2)
                im1 = readFrame(v1);
                im2 = readFrame(v2);
                writeVideo(v,[im1,im2]);
            end
            v.close();

        end
        
        function saveAsVid(filename, im)
            v = VideoWriter(filename);
            v.FrameRate=2;
            v.open();
            for K=1:size(im,3)
                currentSlice = squeeze(im(:,:,K,:));
                writeVideo(v,currentSlice);
            end
            v.close();
        end
        
        function imageStack = readTiff(imname)
            info = imfinfo(imname);
            imageStack = [];
            numberOfImages = length(info);
            for k = 1:numberOfImages
                currentImage = imread(imname, k, 'Info', info);
            end
        end
        
        function newName = getNewFileName(prevFilePath,newDir, ending)
            newName = [newDir prevFilePath(find(prevFilePath=='\',1,'last')+1:end) '.' ending];
        end
        
        function updatedFileStrctCell = translateFileStrctToDifferentPC(fileStrct,oldDir,newDir)
            updatedFileStrctCell = fileStrct;
            
            for ii=1:length(updatedFileStrctCell)
                updatedFileStrctCell{ii}.file = strrep(updatedFileStrctCell{ii}.file,oldDir,newDir);
                updatedFileStrctCell{ii}.seg = strrep(updatedFileStrctCell{ii}.seg,oldDir,newDir);
                for jj=1:length(updatedFileStrctCell{ii}.algoFiles)
                    updatedFileStrctCell{ii}.algoFiles{jj} = strrep(updatedFileStrctCell{ii}.algoFiles{jj},oldDir,newDir);
                end
            end
            
            
        end
        
        function updatedFileStrctCell = saveCroppedFileStrctCell(fileStrctCell, strctName,startPos)
            if ~exist('strctName','var')
                strctName = 'croppedStrctTempName.mat';
            end
            if ~exist('startPos','var')
                startPos = 1;
            end
            if startPos > 1
                updatedFileStrctCell = importdata(strctName);
            else
                updatedFileStrctCell = fileStrctCell;
            end
            updatedFileStrctCell(startPos:end) = [];
            
            for ii=startPos:length(fileStrctCell)
                save(strctName,'updatedFileStrctCell');
                [ii length(fileStrctCell)]
                updatedFileStrctCell{ii} = IO.saveCroppedFileStrctToDir(fileStrctCell{ii});
                
                
                
            end
        end
        
        function copyFilesWithExtensionsRecursively(srcDir,outDir,extension)
            if ~exist(extension,'var')
                extension = 'hd5';
            end
            
            if ~exist(outDir,'dir')
                mkdir(outDir);
            end
            
            files = dir([srcDir '\*.' extension]);
            for file = files'
                inFileName = [srcDir '\' file.name];
                outFileName = [outDir '\' file.name];
                copyfile(inFileName,outFileName);
            end
            
            % Get a list of all files and folders in this folder.
            files = dir(srcDir);
            % Extract only those that are directories.
            subFolders = files([files.isdir]);
            for fold = subFolders'
                if isequal(fold.name,'.') || isequal(fold.name,'..')
                    continue;
                end
                newOutFold = [outDir '\' fold.name];
                newSrcFold = [srcDir '\' fold.name];
                IO.copyFilesWithExtensionsRecursively(newSrcFold,newOutFold,extension);
            end
        end
        
        function [out] = loadFile(fileName)
            if strcmp(fileName(end-2:end),'.gz')
                out = load_untouch_nii_gzip(fileName);
            elseif strcmp(fileName(end-2:end),'hd5')
                out = NiftyToHd5.readHd5(fileName);
            elseif strcmp(fileName(end-2:end),'nii')
                out = load_untouch_nii(fileName);
            else
                error(['input file error: ' fileName]);
            end
        end
        
        function saveFile(obj,fileName)
            if strcmp(fileName(end-2:end),'.gz')
                save_nii_gzip(obj,fileName);
            elseif strcmp(fileName(end-2:end),'hd5')
                error 'not impemented for hs5 yet';
            else
                error(['input file error: ' fileName]);
            end
        end
        
        function updatedFileStrct = saveCroppedFileStrctToDir(fileStrct)
            
            %output
            updatedFileStrct = fileStrct;
            if isempty(fileStrct.seg)
                updatedFileStrct.algoFiles = [];
                return;
            end
            
            %creates new dir
            newImDir = [fileStrct.file(1:find(fileStrct.file=='\',1,'last')) 'CroppedDror\'];
            if ~exist(newImDir, 'dir')
                mkdir(newImDir);
            end
            newSegDir = [fileStrct.seg(1:find(fileStrct.seg=='\',1,'last')) 'CroppedDror\'];
            if ~exist(newSegDir, 'dir')
                mkdir(newSegDir);
            end
            
            
            GT = load_untouch_nii_gzip(fileStrct.seg);
            GT = GT.img;
            im = load_untouch_nii_gzip(fileStrct.file);
            im = im.img;
            
            gap = 25;
            [minY,maxY,minX,maxX,minZ,maxZ] = Utils.extractMinVals(GT,gap);
            
            im2save = im(minY:maxY,minX:maxX,minZ:maxZ);
            newImName = IO.getNewFileName(fileStrct.file,newImDir,'hd5');
            newImName = [newImName(1:end-4) '.' num2str(fileStrct.organ) '.' 'hd5'];
            
            save(newImName, 'im2save', '-v7.3' );
            updatedFileStrct.file = newImName;
            clear im im2save;
            seg2save = GT(minY:maxY,minX:maxX,minZ:maxZ);
            newSegName = IO.getNewFileName(fileStrct.seg,newSegDir,'hd5');
            save(newSegName, 'seg2save', '-v7.3' );
            updatedFileStrct.seg = newSegName;
            clear GT seg2save;
            
            algosToErase = [];
            a
            for t = 1:length(fileStrct.algoFiles)
                [t length(fileStrct.algoFiles)]
                try
                    algo = load_untouch_nii_gzip(fileStrct.algoFiles{t});
                catch err
                    algosToErase = [algosToErase t];
                    continue;
                end
                newAlgoDir = [fileStrct.algoFiles{t}(1:find(fileStrct.algoFiles{t}=='\',1,'last')) 'CroppedDror\'];
                if ~exist(newAlgoDir, 'dir')
                    mkdir(newAlgoDir);
                end
                
                algo = algo.img;
                [minYAlgo,maxYAlgo,minXAlgo,maxXAlgo,minZAlgo,maxZAlgo] = Utils.extractMinVals(algo,gap);
                
                if abs(minYAlgo-minY)>gap || abs(minXAlgo-minX)>gap || abs(minZAlgo-minZ)>gap ...
                        || abs(maxYAlgo-maxY)>gap || abs(maxXAlgo-maxX)>gap || abs(maxZAlgo-maxZ)>gap  ...
                        || isempty(fileStrct.algoFiles{t})
                    algosToErase = [algosToErase t];
                    continue;
                end
                
                algo2save = algo(minY:maxY,minX:maxX,minZ:maxZ);
                newAlgoName = IO.getNewFileName(fileStrct.algoFiles{t},newAlgoDir,'hd5');
                save(newAlgoName, 'algo2save', '-v7.3' );
                updatedFileStrct.algoFiles{t} = newAlgoName;
                
                clear algo algo2save;
            end
            updatedFileStrct.algoFiles(algosToErase) = [];
        end
        function [organFiles] = getRelevantOrganFiles(organNum, directory)
            %fileName = ['1' repmat('0',1,8-length(num2str(num))-1) num2str(num)];
            PATIENT_ID_LENGTH = 8;
            if organNum ==IO.LIVER
                type = '_1_CT_wb';
            else %TODO - type should be parameter
                type = '_1_CTce_ThAb';
            end
            %'*_1_CT_wb'; alternative
            volumeRegex =  ['(.*)' type];
            %volumeRegex =  '*_1_CT_wb'; alternative
            
            
            volumeFiles = regexpdir([directory 'Volumes\'], volumeRegex, false);
            organFiles = cell(length(volumeFiles),1);
            for i=1:length(volumeFiles)
                [i length(volumeFiles)]
                fileStrct.file = volumeFiles{i};
                fileStrct.patientId = fileStrct.file(find(fileStrct.file=='\',1,'last')+1:find(fileStrct.file=='\',1,'last')+1+PATIENT_ID_LENGTH-1);
                fileStrct.organ = organNum;
                segmentationRegex = [fileStrct.patientId '(.*)_' num2str(organNum) '(.*)'];
                algoRegex = [fileStrct.patientId  '(.*)' type '(.*)' num2str(organNum) '(.*)'];
                
                fileStrct.seg = regexpdir([directory 'Segmentations\'], segmentationRegex, false);
                if(length(fileStrct.seg) > 1)
                    segmentationRegex = [fileStrct.patientId '(.*)' type '(.*)' num2str(organNum) '(.*)'];
                    fileStrct.seg = regexpdir([directory 'Segmentations\'], segmentationRegex, false);
                end
                if(length(fileStrct.seg) > 1)
                    error('too much gt for file');
                elseif length(fileStrct.seg) == 1
                    fileStrct.seg = fileStrct.seg{1};
                end
                
                
                fileStrct.algoFiles = regexpdir([directory 'Participants_segments\'], algoRegex, true);
                fileStrct
                organFiles{i} = fileStrct;
            end
            
        end
        
        
        function [out success] = read(directory,num, segNum, algoName, modality)
            PATIENT_ID_LENGTH = 8;
            addpath('NIFTI_analyze_toolbox\');
            addpath('NIFTI_analyze_toolbox_gzip_extension\\');
            fileName = ['1' repmat('0',1,PATIENT_ID_LENGTH-length(num2str(num))-1) num2str(num)];
            fileName = [fileName '_1_CTce_ThAb'];
            
            %reads volume
            if ~exist('segNum','var')
                inputFile = [directory '\Volumes\' fileName '.nii.gz'];
            else
                if ~exist('algoName','var')
                    useSegGT = true;
                else
                    useSegGT = false;
                end
                if ~useSegGT
                    dirName = [directory 'Participants_segments\' algoName '\'];
                else
                    dirName = [directory 'Segmentations\'];
                end
                %reads GT
                segmentationPrefix = [dirName fileName '_' num2str(segNum) '*'];
                inputFile = dir(segmentationPrefix);
                if isempty(inputFile)
                    out = -1;
                    success = false;
                    return;
                end
                
                inputFile =  [dirName inputFile(1).name];
            end
            success = true;
            out = load_untouch_nii_gzip( inputFile);
        end
        
        
        function opt = extractCropOptStrctFromMinMaxVals(volSize,minX,minY,minZ,maxX,maxY,maxZ)
            
            m = volSize(1); n = volSize(2); z = volSize(3);
            opt.cut_from_L = min(int16(m/2-2),minY-1);
            opt.cut_from_R = min(int16(m/2-2), m-maxY);
            opt.cut_from_P = min(int16(n/2-2),minX-1);
            opt.cut_from_A = min(int16(n/2-2),n-maxX);
            opt.cut_from_I = min(int16(z/2-2),minZ-1);
            opt.cut_from_S = min(int16(z/2-2),z-maxZ);
        end
        
        function cropFileStrct(fileStrct)
            gap = 25;
            for ii=1:length(fileStrct)
                [ii length(fileStrct)]
                currStrct = fileStrct{ii};
                
                im = IO.loadFile(currStrct.file);
                
                
                %extracts crop size
                hasSegMask = ~isempty(currStrct.seg);
                
                if hasSegMask
                    seg = IO.loadFile(currStrct.seg);
                    if isstruct(im)
                        [minY,maxY,minX,maxX,minZ,maxZ] = Utils.extractMinVals(seg.img,gap);
                        [m,n,z] = size(seg.img);
                    else
                        [minY,maxY,minX,maxX,minZ,maxZ] = Utils.extractMinVals(seg,gap);
                        [m,n,z] = size(img);
                    end
                elseif isfield(currStrct,'algofiles')
                    minY = inf; minX = inf; minZ = inf;
                    maxY = -inf; maxX = -inf; maxZ = -inf;
                    for jj=1:length(currStrct.algoFiles)
                        try
                            algo = IO.loadFile(currStrct.algoFiles{jj});
                        catch err
                            continue;
                        end
                        if isstruct(algo)
                            algo = algo.img;
                        end
                        [minY2,maxY2,minX2,maxX2,minZ2,maxZ2] = Utils.extractMinVals(algo,gap);
                        [m,n,z] = size(algo);
                        minY = min(minY,minY2); minX = min(minX,minX2); minZ = min(minZ,minZ2);
                        maxY = max(maxY,maxY2); maxX = max(maxX,maxX2); maxZ = max(maxZ,maxZ2);
                    end
                    
                end
                opt = IO.extractCropOptStrctFromMinMaxVals([m,n,z],minX,minY,minZ,maxX,maxY,maxZ);
                %opt.cut_from_L = min(int16(m/2-2),minY-1);
                %opt.cut_from_R = min(int16(m/2-2), m-maxY);
                %opt.cut_from_P = min(int16(n/2-2),minX-1);
                %opt.cut_from_A = min(int16(n/2-2),n-maxX);
                %opt.cut_from_I = min(int16(z/2-2),minZ-1);
                %opt.cut_from_S = min(int16(z/2-2),z-maxZ);
                constClip = false;
                if constClip
                    clipSize = 60;
                    minVal = min(structfun(@(x)min(x(:)),opt));
                    if(minVal<clipSize)
                        error ' cliping is too big';
                    end
                    opt.cut_from_L = clipSize;
                    opt.cut_from_R = clipSize;
                    opt.cut_from_P = clipSize;
                    opt.cut_from_A = clipSize;
                    opt.cut_from_I = clipSize;
                    opt.cut_from_S = clipSize;
                end
                
                if isstruct(im)
                    im = clip_nii(im,opt);
                    if hasSegMask
                        seg = clip_nii(seg,opt);
                    end
                else
                    im = im(minY:maxY,minX:maxX,minZ:maxZ);
                    if hasSegMask
                        seg = seg(minY:maxY,minX:maxX,minZ:maxZ);
                    end
                end
                
                
                algosAndSegs = {currStrct.algoFiles};
                if isfield(currStrct,'manuelSegs')
                    algosAndSegs = [currStrct.algoFiles;currStrct.manuelSegs];
                end
                for jj=1:length(algosAndSegs)
                    try
                        algo = IO.loadFile(algosAndSegs{jj});
                    catch err
                        continue;
                    end
                    if isstruct(im)
                        try
                            algo = clip_nii(algo,opt);
                        catch err
                            warning('algo size is already clipped')
                            size(algo)
                            continue;
                        end
                    else
                        algo = algo(minY:maxY,minX:maxX,minZ:maxZ);
                    end
                    IO.saveFile(algo,algosAndSegs{jj});
                    
                end
                
                
                IO.saveFile(im,currStrct.file);
                if hasSegMask
                    IO.saveFile(seg,currStrct.seg);
                end
                
                clear seg im opt;
                
            end
        end
        function outFileName = copyFileAndGetUpdatedName(fullFileName,folder,appendTopDirName)
            %copyfile(fullFileName,folder);
            
            [oldFolder,name,ext] = fileparts(fullFileName);
            if exist('appendTopDirName','var') && appendTopDirName
                [~,parentFolderName,~] = fileparts(oldFolder);
                outFileName = [folder '\' parentFolderName '_' name ext];
            else
                outFileName = [folder '\' name ext];
            end
            try
                copyfile(fullFileName,outFileName);
            catch
                outFileName = [];
            end
            
        end
        
        function newFileStrct = copyFileStrctToDir(fileStrct, newDir)
            newFileStrct = cell(size(fileStrct));
            for ii=1:length(fileStrct)
                currStrct = fileStrct{ii};
                newFileStrct{ii} = currStrct;
                newFileName = IO.copyFileAndGetUpdatedName(currStrct.file,newDir);
                if ~isempty(currStrct.seg)
                    newSegName = IO.copyFileAndGetUpdatedName(currStrct.seg,newDir);
                end
                newFileStrct{ii}.file = newFileName;
                if ~isempty(currStrct.seg)
                    newFileStrct{ii}.seg = newSegName;
                end
                if isfield(currStrct,'algoFiles')
                    for jj=1:length(currStrct.algoFiles)
                        appendTopDirName = true;
                        newAlgoName = IO.copyFileAndGetUpdatedName(currStrct.algoFiles{jj},newDir, appendTopDirName);
                        newFileStrct{ii}.algoFiles{jj} = newAlgoName;
                    end
                end
                if isfield(currStrct,'manuelSegs')
                    for jj=1:length(currStrct.manuelSegs)
                        appendTopDirName = true;
                        newSegName = IO.copyFileAndGetUpdatedName(currStrct.manuelSegs{jj},newDir, appendTopDirName);
                        newFileStrct{ii}.manuelSegs{jj} = newSegName;
                    end
                end
            end
        end
        
        function  imshow3D( Img, disprange, dim2show )
            %IMSHOW3D displays 3D grayscale images in slice by slice fashion with mouse
            %based slice browsing and window and level adjustment control.
            %
            % Usage:
            % imshow3D ( Image )
            % imshow3D ( Image , [] )
            % imshow3D ( Image , [LOW HIGH] )
            %
            %    Image:      3D image MxNxK (K slices of MxN images)
            %    [LOW HIGH]: display range that controls the display intensity range of
            %                a grayscale image (default: the widest available range)
            %
            % Use the scroll bar or mouse scroll wheel to switch between slices. To
            % adjust window and level values keep the mouse right button pressed and
            % drag the mouse up and down (for level adjustment) or right and left (for
            % window adjustment).
            %
            % "Auto W/L" button adjust the window and level automatically
            %
            % While "Fine Tune" check box is checked the window/level adjustment gets
            % 16 times less sensitive to mouse movement, to make it easier to control
            % display intensity rang.
            %
            % Note: The sensitivity of mouse based window and level adjustment is set
            % based on the user defined display intensity range; the wider the range
            % the more sensitivity to mouse drag.
            %
            %
            %   Example
            %   --------
            %       % Display an image (MRI example)
            %       load mri
            %       Image = squeeze(D);
            %       figure,
            %       imshow3D(Image)
            %
            %       % Display the image, adjust the display range
            %       figure,
            %       imshow3D(Image,[20 100]);
            %
            %   See also IMSHOW.
            
            %
            % - Maysam Shahedi (mshahedi@gmail.com)
            % - Released: 1.0.0   Date: 2013/04/15
            % - Revision: 1.1.0   Date: 2013/04/19
            %
            
            sno = size(Img,3);  % number of slices
            S = round(sno/2);
            
            global InitialCoord;
            global layer;
            
            if ~exist('dim2show','var')
                layer = 3;
            else
                layer = dim2show;
            end
            
            MinV = 0;
            MaxV = max(Img(:));
            LevV = (double( MaxV) + double(MinV)) / 2;
            Win = double(MaxV) - double(MinV);
            WLAdjCoe = (Win + 1)/1024;
            FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients
            
            if isa(Img,'uint8')
                MaxV = uint8(Inf);
                MinV = uint8(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'uint16')
                MaxV = uint16(Inf);
                MinV = uint16(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'uint32')
                MaxV = uint32(Inf);
                MinV = uint32(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'uint64')
                MaxV = uint64(Inf);
                MinV = uint64(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'int8')
                MaxV = int8(Inf);
                MinV = int8(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'int16')
                MaxV = int16(Inf);
                MinV = int16(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'int32')
                MaxV = int32(Inf);
                MinV = int32(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'int64')
                MaxV = int64(Inf);
                MinV = int64(-Inf);
                LevV = (double( MaxV) + double(MinV)) / 2;
                Win = double(MaxV) - double(MinV);
                WLAdjCoe = (Win + 1)/1024;
            elseif isa(Img,'logical')
                MaxV = 0;
                MinV = 1;
                LevV =0.5;
                Win = 1;
                WLAdjCoe = 0.1;
            end
            
            SFntSz = 9;
            LFntSz = 10;
            WFntSz = 10;
            LVFntSz = 9;
            WVFntSz = 9;
            BtnSz = 10;
            ChBxSz = 10;
            
            if (nargin < 2)
                [Rmin Rmax] = WL2R(Win, LevV);
            elseif numel(disprange) == 0
                [Rmin Rmax] = WL2R(Win, LevV);
            else
                LevV = (double(disprange(2)) + double(disprange(1))) / 2;
                Win = double(disprange(2)) - double(disprange(1));
                WLAdjCoe = (Win + 1)/1024;
                [Rmin Rmax] = WL2R(Win, LevV);
            end
            
            %axes('position',[0,0.2,1,0.8]), imshow(Img(:,:,S), [Rmin Rmax])
            if layer==3
                axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])
            elseif layer==2
                axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(:,S,:,:)), [Rmin Rmax])
            else
                axes('position',[0,0.2,1,0.8]), imshow(squeeze(Img(S,:,:,:)), [Rmin Rmax])
            end
            FigPos = get(gcf,'Position');
            S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
            Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
            Wtxt_Pos = [50 20 60 20];
            Wval_Pos = [110 20 60 20];
            Ltxt_Pos = [175 20 45 20];
            Lval_Pos = [220 20 60 20];
            BtnStPnt = uint16(FigPos(3)-250)+1;
            if BtnStPnt < 300
                BtnStPnt = 300;
            end
            Btn_Pos = [BtnStPnt 20 100 20];
            ChBx_Pos = [BtnStPnt+110 20 100 20];
            
            if sno > 1
                shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
                stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
            else
                stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
            end
            ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
            wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
            lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
            wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
            Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
            ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
            
            set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
            set (gcf, 'ButtonDownFcn', @mouseClick);
            set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
            set(gcf,'WindowButtonUpFcn', @mouseRelease)
            set(gcf,'ResizeFcn', @figureResized)
            
            
            % -=< Figure resize callback function >=-
            function figureResized(object, eventdata)
                FigPos = get(gcf,'Position');
                S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
                Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
                BtnStPnt = uint16(FigPos(3)-250)+1;
                if BtnStPnt < 300
                    BtnStPnt = 300;
                end
                Btn_Pos = [BtnStPnt 20 100 20];
                ChBx_Pos = [BtnStPnt+110 20 100 20];
                if sno > 1
                    set(shand,'Position', S_Pos);
                end
                set(stxthand,'Position', Stxt_Pos);
                set(ltxthand,'Position', Ltxt_Pos);
                set(wtxthand,'Position', Wtxt_Pos);
                set(lvalhand,'Position', Lval_Pos);
                set(wvalhand,'Position', Wval_Pos);
                set(Btnhand,'Position', Btn_Pos);
                set(ChBxhand,'Position', ChBx_Pos);
            end
            
            % -=< Slice slider callback function >=-
            function SliceSlider (hObj,event, Img)
                S = round(get(hObj,'Value'));
                %set(get(gca,'children'),'cdata',Img(:,:,S,:))
                if layer==3
                    set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
                elseif layer==2
                    set(get(gca,'children'),'cdata',squeeze(Img(:,S,:,:)))
                else
                    set(get(gca,'children'),'cdata',squeeze(Img(S,:,:,:)))
                end
                
                caxis([Rmin Rmax])
                if sno > 1
                    set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
                else
                    set(stxthand, 'String', '2D image');
                end
            end
            
            % -=< Mouse scroll wheel callback function >=-
            function mouseScroll (object, eventdata)
                UPDN = eventdata.VerticalScrollCount;
                S = S - UPDN;
                if (S < 1)
                    S = 1;
                elseif (S > sno)
                    S = sno;
                end
                if sno > 1
                    set(shand,'Value',S);
                    set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
                else
                    set(stxthand, 'String', '2D image');
                end
                %set(get(gca,'children'),'cdata',Img(:,:,S))
                if layer==3
                    set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
                elseif layer==2
                    set(get(gca,'children'),'cdata',squeeze(Img(:,S,:,:)))
                else
                    set(get(gca,'children'),'cdata',squeeze(Img(S,:,:,:)))
                end
                
            end
            
            % -=< Mouse button released callback function >=-
            function mouseRelease (object,eventdata)
                set(gcf, 'WindowButtonMotionFcn', '')
            end
            
            % -=< Mouse click callback function >=-
            function mouseClick (object, eventdata)
                MouseStat = get(gcbf, 'SelectionType');
                if (MouseStat(1) == 'a')        %   RIGHT CLICK
                    InitialCoord = get(0,'PointerLocation');
                    set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
                end
            end
            
            % -=< Window and level mouse adjustment >=-
            function WinLevAdj(varargin)
                PosDiff = get(0,'PointerLocation') - InitialCoord;
                
                Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
                LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
                if (Win < 1)
                    Win = 1;
                end
                
                [Rmin, Rmax] = WL2R(Win,LevV);
                caxis([Rmin, Rmax])
                set(lvalhand, 'String', sprintf('%6.0f',LevV));
                set(wvalhand, 'String', sprintf('%6.0f',Win));
                InitialCoord = get(0,'PointerLocation');
            end
            
            % -=< Window and level text adjustment >=-
            function WinLevChanged(varargin)
                
                LevV = str2double(get(lvalhand, 'string'));
                Win = str2double(get(wvalhand, 'string'));
                if (Win < 1)
                    Win = 1;
                end
                
                [Rmin, Rmax] = WL2R(Win,LevV);
                caxis([Rmin, Rmax])
            end
            
            % -=< Window and level to range conversion >=-
            function [Rmn Rmx] = WL2R(W,L)
                Rmn = L - (W/2);
                Rmx = L + (W/2);
                if (Rmn >= Rmx)
                    Rmx = Rmn + 1;
                end
            end
            
            % -=< Window and level auto adjustment callback function >=-
            function AutoAdjust(object,eventdata)
                Win = double(max(Img(:))-min(Img(:)));
                Win (Win < 1) = 1;
                LevV = double(min(Img(:)) + (Win/2));
                [Rmin, Rmax] = WL2R(Win,LevV);
                caxis([Rmin, Rmax])
                set(lvalhand, 'String', sprintf('%6.0f',LevV));
                set(wvalhand, 'String', sprintf('%6.0f',Win));
            end
            
        end
        % -=< Maysam Shahedi (mshahedi@gmail.com), April 19, 2013>=-
    end
    
end

