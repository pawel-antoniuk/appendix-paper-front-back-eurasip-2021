% Paweł Antoniuk 2021
% Bialystok University of Technology

classdef H5MatrixDatastore < matlab.io.Datastore & ...
        matlab.io.datastore.Shuffleable & ...
        matlab.io.datastore.Partitionable
    properties (Access = private)
        CurrentIndex double
        Filename string        
        DataLocation string
        LabelsLocation string
        MatrixSize                
    end
    
    properties (Access = public)        
        IndexValues
        Labels
        NumMatrices
    end
    
    methods
        function ds = H5MatrixDatastore(filename, dataLocation, labelsLocation)
            ds.CurrentIndex = 1;
            ds.Filename = filename;
            ds.DataLocation = dataLocation;
            ds.LabelsLocation = labelsLocation;
            hsize = h5info(ds.Filename, ds.DataLocation).Dataspace.Size;
            ds.MatrixSize = hsize(1: end-1);
            ds.NumMatrices = hsize(end);
            ds.IndexValues = randperm(ds.NumMatrices);
            ds.Labels = categorical(h5read(ds.Filename, ds.LabelsLocation));
%             ds.Labels = h5read(ds.Filename, ds.LabelsLocation);
            
            reset(ds);
        end
        
        function tf = hasdata(ds)
            tf = ds.CurrentIndex <= ds.NumMatrices;
        end
        
        function [data, info] = read(ds)
            if ~hasdata(ds)
                error('The H5 file has no more data');
            end
            
            indx = ds.IndexValues(ds.CurrentIndex);
            
            dataBegin = [ones(1, numel(ds.MatrixSize)), indx];
            dataCount = [ds.MatrixSize, 1];
            data{1} = h5read(ds.Filename, ds.DataLocation, dataBegin, dataCount);
            data{2} = ds.Labels(indx);
            
            info.Size = size(data);
            info.FileName = [ds.Filename, '/', ds.DataLocation, indx];
            
            ds.CurrentIndex = ds.CurrentIndex + 1;
        end
        
        function reset(ds)
            ds.CurrentIndex = 1;
        end
        
        function dsNew = shuffle(ds)
            dsNew = copy(ds);
            dsNew.IndexValues = dsNew.IndexValues(randperm(dsNew.NumMatrices));
        end
        
        function subds = partition(ds, n, index)
            partitionSize =  ceil(ds.NumMatrices / n);
            indxBegin = (index - 1 ) * partitionSize + 1;
            indxEnd = min(index * partitionSize, ds.NumMatrices);
            
            subds = copy(ds);

            subds.IndexValues = subds.IndexValues(indxBegin:indxEnd);
            subds.NumMatrices = length(subds.IndexValues);
            reset(subds);
        end
        
        function [subds1, subds2] = splitDataset(ds, p)
            splitPoint = round(numel(ds.IndexValues) * p);
            subds1 = copy(ds);
            subds1.IndexValues = ds.IndexValues(1:splitPoint);
            subds1.NumMatrices = numel(subds1.IndexValues);
            
            subds2 = copy(ds);
            subds2.IndexValues = ds.IndexValues(splitPoint+1:end);
            subds2.NumMatrices = numel(subds2.IndexValues);
        end
        
        function [filteredDs, restDs] = filterDs(ds, filterFunc)
            filteredDs = copy(ds);
            filteredDs.IndexValues = [];
            restDs = copy(ds);
            restDs.IndexValues = [];
            
            for ii = ds.IndexValues
                if filterFunc(ii)
                    filteredDs.IndexValues(end + 1) = ii;
                else
                    restDs.IndexValues(end + 1) = ii;
                end
            end
            
            filteredDs.NumMatrices = numel(filteredDs.IndexValues);
            restDs.NumMatrices = numel(restDs.IndexValues);
        end
    end
    
    methods (Access = protected)
        function n = maxpartitions(ds)
            n = ds.NumMatrices;
        end
    end
    
    methods (Hidden = true)
        function frac = progress(ds)
            if hasdata(ds)
                frac = (ds.CurrentIndex - 1) /  ds.NumMatrices;
            else
                frac = 1;
            end
        end
    end
end






