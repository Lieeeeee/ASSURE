    classdef SUREExperimentCls
    %SUREEXPERIMENTCLS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        livers = {'FU1','FU2','FU4','FU5','FU7'};
        lungs = {'case3','case4','case5','case10','case11'};
        kidneys = {'10000100','10000104','10000105','10000109','10000111','10000112'};
        brains = {'caseb2','caseb5'};
        
        longDirNames = {'largeAmount'};
        novicesDirs = {'largeAmount','mediumAmount','novices'};
        shortDirNames = {'largeAmount','mediumAmount'};
        extraShortDirNames = {'largeAmount','mediumAmount','smallAmount'};
        
        slicesExtraShortLivers = {[15,20,25,31,36,41], [],[], [35,39,43,94,95,96], []};
        slicesExtraShortLungs = { [],[],[120,123,126,130], [30,34], [123,129,130]};
        slicesExtraShortKydneys = {  [217,226],[],[],[192,203,208,215],[197,202],[190,203,211,217]}
        
        slicesShortLivers = {[15,20,25,31,36,41], [72,74,77,79], [31,33,35,38 ], [31,35,39,43,94,95,96,99], [114,118,122,126,130]};
        slicesShortLungs = {[55, 60, 65, 72, 76], [44,48,52,56, 60], [120,123,126,130], [30,32,34,35], [110,115,119,123,129,130]};
        slicesShortKydneys = {  [190,200,210,217,226], [207,221,235,249], [200,210,220,230],[192,203,208,215,221],[197,202,212,222,227],[190,195,203,211,217]}
        slicesShortBrains = {  [40,50,60,65,70,75], [23,27,31,33, 37]};
        
        slicesShort = [SUREExperimentCls.slicesShortLivers(:);SUREExperimentCls.slicesShortLungs(:);  SUREExperimentCls.slicesShortKydneys(:)];
        
        sliceLongLivers = {10:44, 69:81, 31:38, [27:46,90:99],111:136}; %35+13+8+20+10+26
        sliceLongBrains = {33:76,19:42};
        %sliceLongLivers = {7:47, 67:81, 30:39, [26:46,90:101],109:136};
        %sliceLongLungs = {41:83,30:66, 113:139,27:36,103:140};
        sliceLongLungs = {43:81,30:66, 114:139,28:36,105:140}; %43+37+27+10+38
        sliceLongKydneys = { ...
            [170,175,180,185,190,195,200,205,210,215,217,220,226], ... %10000100
            [207,214,221,228,235,242,249,256,263,270,277], ... %10000104
            [190,200,210,220,230,240], ...10000105
            [192,195,198,201,203,208,210,213,215,218,221], ... 10000109
            [197,202,207,212,217,222,227,232], ... %10000111
            [186,190,192,195,198,201,203,207,211,213,217,219,222]}; %10000112
        
        slicesLong = [SUREExperimentCls.sliceLongLivers(:);SUREExperimentCls.sliceLongLungs(:);  SUREExperimentCls.sliceLongKydneys(:)];
        %slicesLong = { ...
        %    7:47, 67:81, 30:39, [26:46,90:101],109:136, ... %livers
        %    41:83,30:66, 113:139,27:36,103:140, ... %lungs
        %    [170,175,180,185,190,195,200,205,210,215,217,220,226], ... %10000100
        %    [207,214,221,228,235,242,249,256,263,270,277], ... %10000104
        %    [190,200,210,220,230,240], ...10000105
        %    [192,195,198,201,203,208,210,213,215,218,221], ... 10000109
        %    [197,202,207,212,217,222,227,232], ... %10000111
        %    [186,190,192,195,198,201,203,207,211,213,217,219,222]}; %10000112
        
        windowingLivers = {[-130,220],[-130,220],[-80,80],[-130,220],[-80,150]};
        
        windowingLungs = {[-1024,-80],[-1024,60],[-1024,415],[-11262,2000],[-1024,180]};
        
        windowingKidneys = {[-130,220],[-130,220],[-200,250],[-200,300],[-200,200],[-150,200]};
        
        windowingBrains = {[995,1124],[995,1150]};
        
        windowing = {[-130,220],[-130,220],[-80,80],[-130,220],[-80,150],[-1024,-80],[-1024,60],[-1024,415],[-11262,2000],[-1024,180], ...
            [-130,220],[-130,220],[-200,250],[-200,300],[-200,200],[-150,200],[995,1124],[995,1150]};
    end
    
    methods (Static)
        
        
        function [ sen] = getRadiologistSeniority(str)
            if length(findstr('as',str))>0
                sen = 2;
            elseif length(findstr('dh',str))>0
                sen = 2;
            elseif length(findstr('bg',str))>0 %ben gurion -zahi annotations
                sen = 2;
            elseif length(findstr('bd',str))>0
                sen = 4;
            elseif length(findstr('nl',str))>0
                sen = 4;
            elseif length(findstr('NC',str))>0
                sen = 4;
            elseif length(findstr('ns',str))>0  || length(findstr('ns2',str))>0
                sen = 4;
            elseif length(findstr('ma',str))>0
                sen = 3;
            elseif length(findstr('kh',str))>0
                sen = 3;
            elseif length(findstr('sb',str))>0 || length(findstr('sb2',str))>0
                sen = 3;
            elseif length(findstr('js',str))>0
                sen = 4;
            elseif length(findstr('ok',str))>0 %Onur alp karakasalan
                sen = 3;
            elseif length(findstr('nl1',str))>0
                sen = 4;
            elseif length(findstr('ok',str))>0
                sen = 3;
            elseif length(findstr('gt',str))>0
                sen = 4;
            elseif length(findstr('az',str))>0 || length(findstr('ch',str))>0
                sen = 1;
            else
                error 'didnt find radiologists'
            end
        end
        function [name, sen] = getRadiologistName(fileName)
            if length(findstr('SOTO',fileName))>0
                name = 'as';
            elseif length(findstr('_dh',fileName))>0
                name = 'dh';
            elseif length(findstr('_BD',fileName))>0
                name = 'bd';
                
            elseif length(findstr('_NC',fileName))>0
                name = 'NC';
            elseif length(findstr('natasha2',fileName))>0
                name = 'ns2';
            elseif length(findstr('natasha',fileName))>0
                name = 'ns';
                
            elseif length(findstr('AWAD',fileName))>0
                name = 'ma';
            elseif length(findstr('khaled',fileName))>0
                name = 'kh';
            elseif length(findstr('sasha2',fileName))>0
                name = 'sb2';
            elseif length(findstr('sasha',fileName))>0
                name = 'sb';
                
            elseif length(findstr('Sosna',fileName))>0 || length(findstr('sosna',fileName))>0
                name = 'js';
            elseif length(findstr('onur',fileName))>0 %Onur alp karakasalan
                name = 'ok';
            elseif length(findstr('naama1',fileName))>0
                name = 'nl1';
            elseif length(findstr('_nlc',fileName))>0 || length(findstr('naama',fileName))>0
                name = 'nl';
            elseif length(findstr('gtlungs',fileName))>0
                name = 'ok';
            elseif length(findstr('bench',fileName))>0
                name = 'gt';
            elseif length(findstr('bgusemi',fileName))>0 %ben gurion -zahi annotations
                name = 'bg1';
            elseif length(findstr('bgu',fileName))>0 %ben gurion -zahi annotations
                name = 'bg';
            elseif length(findstr('clara',fileName))>0
                name = 'ch';
            elseif length(findstr('az',fileName))>0
                name = 'az';
                
            else
                error 'didnt find radiologists'
            end
        end
    end
    
end

