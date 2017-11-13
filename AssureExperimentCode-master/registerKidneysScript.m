%cropps files
%filesStrctCell = importdata('matStructs\kidneyLeftST.mat');
%newDir = 'C:\Users\drorcohe\Desktop\thesisWork\leftKidneyFiles\ST';
%newFileStrct = IO.copyFileStrctToDir(filesStrctCell(1:6), newDir);
%save('matStructs\kidneyLeftSTMovedPartial.mat','newFileStrct');
newFileStrct = importdata('matStructs\kidneyLeftSTMovedPartial.mat');
%newFileStrct = IO.translateFileStrctToDifferentPC(newFileStrct,'C:\Users\drorcohe\Desktop\thesisWork\','C:\Users\drorcohe\Desktop\thesisProj\segmentationevaluation\');

%IO.cropFileStrct(newFileStrct);

%generates Atlas
gtStrctCell = importdata('matStructs\kidneyLeftGTMoved.mat');
%gtStrctCell = IO.translateFileStrctToDifferentPC(gtStrctCell,'C:\Users\drorcohe\Desktop\thesisWork\','C:\Users\drorcohe\Desktop\SURECodeBase\SUREModel\');

%gtStrctCell = IO.translateFileStrctToDifferentPC(gtStrctCell,'C:\Users\drorcohe\Desktop\thesisWork\','C:\Users\drorcohe\Desktop\thesisProj\segmentationevaluation\');
opened = [2,4,7,10,13,14,17,19,20];
extremelyOpen = [12];
closed = [1,5,6, 8, 9,11,15,16,18];
invalid = [3];
atlas = ShapeAtlas(gtStrctCell(opened));
atlas = atlas.initAtlas();

%calc shape priot by Atlas
imNii = IO.loadFile(newFileStrct{1}.file);
segNii = IO.loadFile(newFileStrct{1}.seg);
algoSegNii = IO.loadFile(newFileStrct{1}.algoFiles{1});
[res, shapePriorMask] = Prior.shapePriorByAtlas(imNii,algoSegNii,atlas);
dispRes = Utils.displayCertaintyUncertainty2_3D(mat2gray(imNii.img),res,Utils.getBoundries(algoSegNii.img));
dispRes2 = Utils.displaySegmentation(dispRes,Utils.getBoundries(shapePriorMask),[0,0,1]);

bbox1 = Utils.getBoundingBox(segNii.img);
bbox2 = Utils.getBoundingBox(algoSegNii.img);

imRefNii = IO.loadFile(newFileStrct{3}.file);
segRefNii = IO.loadFile(newFileStrct{3}.seg);
segRef2Nii = IO.loadFile(newFileStrct{3}.algoFiles{1});
[out, seg] = NiftyReg.registerTwoNiis(imNii,imRefNii,segNii,segRefNii);
[out2, seg2] = NiftyReg.registerTwoNiis(imNii,imRefNii,segNii,segRef2Nii);


%%
%registers images in file strct to the first file.


%filesStrctCell = importdata('matStructs\kidneyLeftGT.mat');
%newDir = 'C:\Users\drorcohe\Desktop\thesisWork\leftKidneyFiles';
%newFileStrct = IO.copyFileStrctToDir(filesStrctCell, newDir);
filesStrctCell = importdata('matStructs\kidneyLeftGTMovedCropped.mat');
%filesStrctCell = IO.translateFileStrctToDifferentPC(filesStrctCell,'C:\Users\drorcohe\Desktop\thesisWork\','C:\Users\drorcohe\Desktop\thesisProj\segmentationevaluation\');

opened = [2,4,7,10,13,14,17,19,20];
extremelyOpen = [12];
closed = [1,5,6, 8, 9,11,15,16,18];
invalid = [3];
%IO.cropFileStrct(newFileStrct);
[filesStrctCellAlignedOpen, res]= NiftyReg.registerFileStrct(filesStrctCell(opened),'leftKidneyFiles\aligned');
[filesStrctCellAlignedClose, res]= NiftyReg.registerFileStrct(filesStrctCell(closed),'leftKidneyFiles\aligned');

filesStrctCellAligned = filesStrctCellAlignedClose;
niiImages = cell(size(filesStrctCellAligned));
niiMasks = cell(size(filesStrctCellAligned));
segHist = [];
for ii=1:length(filesStrctCellAligned)
    im = IO.loadFile(filesStrctCellAligned{ii}.file);
    niiImages{ii} = im.img;
    seg = IO.loadFile(filesStrctCellAligned{ii}.seg);
    seg.img = seg.img > 0;
    niiMasks{ii} = seg.img;
    if isempty(segHist)
       segHist = seg.img;
    elseif res(ii)==1
        segHist = segHist + seg.img;
    end
end

NiftyReg.registerFileStrct(filesStrctCell{open},'C:\Users\drorcohe\Desktop\thesisProj\segmentationevaluation\leftKidneyFiles\aligned');



openStrct = filesStrctCell{open};

registerImages(openStrct);