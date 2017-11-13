classdef NiftyReg
    %NIFTYREG Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        %APPINE_T_PATH = 'C:\Users\drorcohe\Desktop\build\reg-apps\Release\reg_aladin.exe';
        %DEFORMANLE_T_PATH = 'C:\Users\drorcohe\Desktop\build\reg-apps\Release\reg_f3d.exe';
        %JACOBIAN_T_PATH = 'C:\Users\drorcohe\Desktop\build\reg-apps\Release\reg_resample.exe';
        APPINE_T_PATH = 'C:\Users\drorcohe\Desktop\nfdd\build\reg-apps\Release\reg_aladin.exe';
        DEFORMANLE_T_PATH = 'C:\Users\drorcohe\Desktop\nfdd\build\reg-apps\Release\reg_f3d.exe';
        JACOBIAN_T_PATH = 'C:\Users\drorcohe\Desktop\nfdd\build\reg-apps\Release\reg_resample.exe';
        
    end
    properties
        
        
    end
    
    methods (Static)
        
        function [imOut,segOut,extraSegToTransform] = registerTwoNiis(im,imRef,seg,segRef,segToTransform)
            warning('this function uses unnecessarry read and write. it is for debug')
            bboxSeg = seg;
            bboxSeg.img = Utils.getBoundingBox(bboxSeg.img>0);
            bboxSegRef = segRef;
            bboxSegRef.img = Utils.getBoundingBox(bboxSegRef.img>0);
            
            TEMP_IM_NAME = 'temp\imNii.nii.gz';
            TEMP_SEG_NAME = 'temp\segNii.nii.gz';
            TEMP_IM_REF_NAME = 'temp\imRefNii.nii.gz';
            TEMP_SEG_REF_NAME = 'temp\segRefNii.nii.gz';
            IM_OUT_NAME = 'temp\outImNii.nii.gz';
            SEG_OUT_NAME = 'temp\outSegNii.nii.gz';
            im.untouch = 0;
            seg.untouch = 0;
            imRef.untouch = 0;
            segRef.untouch = 0;
            bboxSeg.untouch = 0;
            bboxSegRef.untouch = 0;
            save_nii_gzip(im,TEMP_IM_NAME);
            %save_nii_gzip(seg,TEMP_SEG_NAME);
            save_nii_gzip(bboxSeg,TEMP_SEG_NAME);
            save_nii_gzip(imRef,TEMP_IM_REF_NAME);
            %save_nii_gzip(segRef,TEMP_SEG_REF_NAME);
            save_nii_gzip(bboxSegRef,TEMP_SEG_REF_NAME);
            
            NiftyReg.registerNiftysAndMasks(TEMP_IM_REF_NAME,TEMP_IM_NAME,TEMP_SEG_REF_NAME,TEMP_SEG_NAME,IM_OUT_NAME,SEG_OUT_NAME);
            delete TEMP_IM_NAME;
            delete TEMP_SEG_NAME;
            delete TEMP_IM_REF_NAME;
            delete TEMP_SEG_REF_NAME;
            
            segOut = IO.loadFile(SEG_OUT_NAME);
            imOut = IO.loadFile(IM_OUT_NAME);
            delete SEG_OUT_NAME;
            delete IM_OUT_NAME;
        end
        
        
        function [fileStrctOut resVec] = registerFileStrct(fileStrct,newFolder)
            fileStrctOut = fileStrct;
            firstNii = fileStrct{1}.file;
            firstNiiSeg = fileStrct{1}.seg;
            
            resVec = zeros(length(fileStrct),1);
            resVec(1) = 1;
            for ii=2:length(fileStrct)
                [ii length(fileStrct)]
                currFileStrct = fileStrct{ii};
                fileStrctOut{ii} = currFileStrct;
                [pathstr,name,ext] = fileparts(currFileStrct.file);
                newNiiName = [newFolder '\' name ext];
                [pathstr,name,ext] = fileparts(currFileStrct.seg);
                newSegName = [newFolder '\' name ext];
                if(ii==1)
                    copyfile(currFileStrct.file,newNiiName);
                    copyfile(currFileStrct.seg,newSegName);
                    
                else
                    res = NiftyReg.registerNiftysAndMasks(firstNii,currFileStrct.file,firstNiiSeg,currFileStrct.seg, newNiiName, newSegName);
                    resVec(ii) = res;
                end
                fileStrctOut{ii}.file = newNiiName;
                fileStrctOut{ii}.seg = newSegName;
                
            end
        end
        
        function res = registerNiftysFromTransformation(niiRef,niiFlo,niiOut, cppFileName)
            masksInputStr = [' -ref ' niiRef ' -flo ' niiFlo];
            system([NiftyReg.JACOBIAN_T_PATH  masksInputStr ' -cpp ' cppFileName  ' -res ' niiOut]);
        end
        
        function res = registerNiftysAndMasks(niiRef,niiFlo,segRef,segFlo, out, outSeg, tranformationFileName)
            if ~exist(NiftyReg.APPINE_T_PATH)==2
                error('NIFTY Reg not found');
            end
            if exist('tranformationFileName','var')
                cppFileName = tranformationFileName;
            else
                cppFileName = 'cpp.nii';
            end
            
            niiInputStr = [' -ref ' niiRef ' -flo ' niiFlo];
            maskInputStr = [ ' -rmask ' segRef ' -fmask ' segFlo];
            masksInputStr = [' -ref ' segRef ' -flo ' segFlo];
            %perform affine transofmation on image
            %clc
            %fprintf('affine reg stage \n');
            system([NiftyReg.APPINE_T_PATH niiInputStr maskInputStr ' -res outAffRes.nii' ]);
            %perform affine transofmation on mask
            %clc
            %fprintf('deform reg stage \n');
            system([NiftyReg.DEFORMANLE_T_PATH niiInputStr maskInputStr ' -aff outputAffine.txt -res outDeform.nii -cpp ' cppFileName]);
            %apply transformation (not required)
            system([NiftyReg.JACOBIAN_T_PATH  niiInputStr ' -cpp ' cppFileName  ' -res ' out]);
            %apply transformation on masks
            system([NiftyReg.JACOBIAN_T_PATH  masksInputStr ' -cpp ' cppFileName  ' -res ' outSeg]);
            
            %clc
            fprintf('clear env \n');
            if ~exist('tranformationFileName','var')
               delete('cpp.nii');
            end
            delete('outputAffine.txt');
            delete('outDeform.nii');
            delete('outAffRes.nii');
            
            %clc
            %fprintf('display results\n');
            regS = load_nii_gzip(niiRef);
            floS = load_nii_gzip(niiFlo);
            outRes = load_nii_gzip(out);
            
            %figure,IO.imshow3D(regS.img);
            %%   figure,IO.imshow3D(floS.img);
            %figure,IO.imshow3D(outRes.img);
            close all;
            allIm = cat(2,mat2gray(outRes.img),mat2gray(regS.img));
            flo2show = NiftyReg.padIm(floS.img,size(allIm));
            allIm = cat(2,mat2gray(flo2show),allIm);
            %figure,IO.imshow3D(floS.img);
            figure,IO.imshow3D(allIm);
            %in = input('is it ok? press 0 if not\n');
            %res = in ~= 0;
            %if res==0
            %    a = 3;
            %end
            
        end
        
        function res = padIm(im,newSize)
            res = im;
            if size(res,1)>newSize(1)
                res = res(1:newSize(1),:,:);
            else
                res = [res;zeros(newSize(1)-size(res,1),size(res,2),size(res,3))];
            end
            if size(res,3)>newSize(3)
                res = res(:,:,1:newSize(3));
            else
                res = cat(3,res,zeros(size(res,1),size(res,2),newSize(3)-size(res,3)));
            end
        end
    end
    
    
    
end

