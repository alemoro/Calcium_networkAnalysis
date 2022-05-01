classdef networkActivityApp < matlab.apps.AppBase
    %NETWORKACTIVITYAPP is used to analyze Ca2+ imaging movies in a
    %network with one extra image as information on cells (for now DIC)
    %   Detailed explanation goes here
    
    % app components
    properties (Access = public)
        Figure
        FileMenu
        FileMenuOpen
        FileMenuLoad
        FileMenuSave
        FileMenuExport
        FileMenuLabelCondition
        DicMenu
        DicMenuPlaceRoi
        DicMenuRemoveRoi
        DicMenuSelect
        AnalysisMenu
        AnalysisMenuDetect
        AnalysisMenuQuantify
        AnalysisMenuReCalculate
        PlotMenu
        PlotMenuRaster
        PlotMenuRasterAll
        PlotMenuRasterSpike
        PlotMenuRasterTrace
        PlotMenuOverview
        PlotMenuTimeFrequency
        OptionMenu
        OptionMenuSettings
        OptionMenuDebug
        AxesDIC
        bcDIC
        AxesStack
        manualReg
        AxesPlot
        ToggleViewStack
        SliderImageStack
        ListImagesLabel
        ListImages
        ButtonPlot
        RadioAllMean
        RadioSingleTrace
        CellNumberLabel
        CellNumberText
        CellNumberPlus
        CellNumberMinor
        TogglePeakAdd
        TogglePeakRemove
        ToggleFixAxis
    end
    
    % File storage properties
    properties (Access = private)
        imgDatastore % Collection of all Images
        stackDs % datastore collection for the tall array that keeps all the stack images
        dicT % Table containing the information from the DIC image
        imgT % Table containing the information of the peaks
        fullT % table with the imgT plus the additional informations on the peaks
        cellsMask % B&W image as mask of the cells ROIs
        patchMask % container for the cells mask
        curDIC % The DIC image that is display
        hDIC % Handle to DIC image
        curStak % The stack that is selected
        hStack % Handle to movie image
        minMax % Minimum and Maximum value for DIC and movie
        curCell % The cell that is selected
        movieData % Temporary storage of the timelapse
        curSlice % The timelapse slice that is visible
        curTime % a line for the current position on the plot
        options % Store the options
        YMinMax % Store the minimum and maximum value of the Y axis in this FOV
    end
    
    % Interaction methods
    methods (Access = private)
        function togglePointer(app)
            if strcmp(app.Figure.Pointer, 'arrow')
                app.Figure.Pointer = 'watch';
            else
                app.Figure.Pointer = 'arrow';
            end
            drawnow();
        end
              
        % Populate the Image list
        function populateImageList(app)
            if isempty(app.dicT) || sum(contains(app.dicT.CellID, app.curDIC)) ~= 1
                dicIDs = regexp(app.curDIC(1:end-4), '_', 'split');
                dicID = [dicIDs{1} '_' dicIDs{2}];
            else
                dicID = app.dicT.ExperimentID{contains(app.dicT.CellID, app.curDIC)};
            end
            app.ListImages.String = app.imgT{strcmp(app.imgT.ExperimentID, dicID), 'CellID'};
            if isempty(app.ListImages.String)
                return
            end
            % if there is data plot the data
            [~, curStack] = fileparts(app.imgT.Filename{strcmp(app.imgT.CellID,app.ListImages.String{1})});
            app.curStak = app.ListImages.String{1};
            if ~isempty(app.imgT)
                if ~isempty(cell2mat((app.imgT{strcmp(app.imgT.CellID, app.curStak),'RawIntensity'})))
                    updatePlot(app);
                else
                    cla(app.AxesPlot)
                end
            end
        end
        
        % Create a B&W image to store the cell ROIs
        function getCellMask(app, whatDIC, sizeImage)
            r = app.dicT.RoiSet{whatDIC}(:,1);
            c = app.dicT.RoiSet{whatDIC}(:,2);
            % not for now since it merges connected regions
            %             % create a mask to identify the neurons
            %             dicMask = false(sizeImage);
            %             for ij = 1:numel(r)
            %                 dicMask(r(ij)-5:r(ij)+5, c(ij)-5:c(ij)+5) = true;
            %             end
            %             app.cellsMask = bwboundaries(dicMask);
            app.cellsMask = cell(numel(r), 2);
            for ij = 1:numel(r)
                app.cellsMask{ij,1} = [r(ij)-app.options.RoiSize r(ij)+app.options.RoiSize r(ij)+app.options.RoiSize r(ij)-app.options.RoiSize];
                app.cellsMask{ij,2} = [c(ij)-app.options.RoiSize c(ij)-app.options.RoiSize c(ij)+app.options.RoiSize c(ij)+app.options.RoiSize];
            end
        end
        
        % Load the DIC Image on the axes
        function updateDIC(app, bTable)
            app.AxesDIC.Visible = 'on';
            app.bcDIC.Visible = 'on';
            % Check if DIC is already present
            if bTable
                if isempty(app.curDIC)
                    app.curDIC = app.dicT.CellID{1};
                end
                whatDIC = contains(app.dicT.CellID, app.curDIC);
            else
                whatDIC = 0;
            end
            cellRoi = [];
            if sum(whatDIC) == 1
                dicName = regexprep(app.curDIC, '_', ' ');
                dicImage2 = app.dicT.CorImage{whatDIC};
                app.minMax.DIC = stretchlim(dicImage2);
                if ~isempty(app.dicT.RoiSet{whatDIC})
                    getCellMask(app, whatDIC, size(dicImage2))
                    cellRoi = app.cellsMask;
                end
            else
                dicName = regexprep(app.curDIC, '_', ' ');
                dicFile = app.dicT.Filename(contains(app.dicT.Filename, app.curDIC));
                dicImage = imread(dicFile{:});
                % Correct DIC image for bad illumination
                dicImage2 = imcomplement(dicImage);
                filterSE = strel('square', 10);
                dicBackground = imerode(dicImage2, filterSE);
                dicBackground = imdilate(dicBackground, filterSE);
                dicBackground = imgaussfilt(dicBackground, 10);
                dicBackground = dicBackground - min(min(dicBackground));
                dicImage2 = imcomplement(dicImage2 - dicBackground);
                % For now detect manually the neurons
                app.minMax.DIC = stretchlim(dicImage2);
            end
            if app.imgT.ImgProperties(1) < size(dicImage2,1)
                dicImage2 = imresize(dicImage2, [app.imgT.ImgProperties(1) app.imgT.ImgProperties(2)]);
            end
            app.hDIC = imshow(imadjust(dicImage2, app.minMax.DIC), 'Parent', app.AxesDIC);
            app.AxesDIC.XLim = [0 size(dicImage2, 1)];
            app.AxesDIC.YLim = [0 size(dicImage2, 2)];
            app.AxesDIC.Title.String = dicName;
            % superimpose the ROI if any
            if ~isempty(cellRoi)
                showMaskDIC(app, cellRoi);
            else
                app.cellsMask = [];
            end
            populateImageList(app)
            app.AxesDIC.Toolbar.Visible = 'off';
        end
        
        function bcDICpressed(app, event)
            whatDIC = contains(app.dicT.CellID, app.curDIC);
            dicImage2 = app.dicT.CorImage{whatDIC};
            stackImg = app.hStack.CData;
            if size(stackImg,1) < size(dicImage2,1)
                dicImage2 = imresize(dicImage2, 0.5);
            end
            adjustBrightness(app, dicImage2, stackImg);
        end
        
        % Load and display the movie, not stored localy, so prone to error but faster for the GUI
        function createStackMovie(app)
            warning('off', 'all');
            % Get the current file
            imgFile = app.imgT.Filename{contains(app.imgT.Filename, app.curStak)};
            imgInfo = app.imgT{contains(app.imgT.CellID, app.curStak),'ImgProperties'};
            imgWidth = imgInfo(1);
            imgHeight = imgInfo(2);
            imgNumber = imgInfo(3);
            imgStack = zeros(imgHeight, imgWidth, imgNumber, 'uint16');
            %startPoint = cell2mat(app.imgT{contains(app.imgT.CellID, app.curStak),'ImgByteStrip'});
            imgStackID = fopen(imgFile, 'r');
            tstack = Tiff(imgFile);
            imgStack(:,:,1) = tstack.read();
            %hWait = waitbar(0, 'Loading image...');
            hWait = waitbar(0, 'Loading image without progress. It''s faster!');
            imgStack(:,:,1) = tstack.read();
            for n = 2:imgNumber
                % Check for clicked Cancel button
                %waitbar(n/imgNumber, hWait, sprintf('Loading image "%s" (%0.2f%%)', regexprep(app.curStak, '_', ' '), n/imgNumber*100));
                nextDirectory(tstack);
                imgStack(:,:,n) = tstack.read();
                %fseek(imgStackID, startPoint(n), 'bof');
                %imgStack(:,:,n) = fread(imgStackID, [imgWidth imgHeight], 'uint16=>uint16');
            end
            fclose(imgStackID);
            delete(hWait)
            % If present add the ROIs
            
            % Adjust the orientation of the image
            %imgStack = flip(imgStack, 1);
            %imgStack = imrotate(imgStack, -90);
            % Data might need registration
            [~, fileName] = fileparts(imgFile);
            fileFltr = find(contains(app.imgT.CellID, fileName));
            if ~contains(fileName, app.options.Reference) && app.options.Registration
                if any(strcmp('RegEst', app.imgT.Properties.VariableNames))
                    if isempty(cell2mat(app.imgT{fileFltr, 'RegEst'}))
                        tyrFile = app.imgT{find(contains(app.imgT{1:fileFltr,'CellID'}, app.options.Reference),1,'last'),'Filename'};
                        imgTyr = imread(cell2mat(tyrFile));
                        imgOut = imref2d(size(imgTyr));
                        imgTif = imgStack(:,:,1);
                        tformEstimate = imregcorr(imgTif, imgTyr, 'translation');
                        app.imgT{fileFltr,'RegEst'} = {tformEstimate.T};
                    else
                        tformEstimate = affine2d(cell2mat(app.imgT{fileFltr,'RegEst'}));
                        imgOut = imref2d(size(imgStack(:,:,1)));
                    end
                end
                imgStack = imwarp(imgStack, tformEstimate, 'OutputView', imgOut);
            end
            app.movieData = imgStack;
            app.SliderImageStack.Visible = 'on';
            app.SliderImageStack.Value = 1;
            app.SliderImageStack.SliderStep = [1/(imgNumber-1), 1/(imgNumber-1)];
            app.SliderImageStack.Min = 1;
            app.SliderImageStack.Max = imgNumber;
            s = imshow(app.movieData(:,:,1), [min(min(min(app.movieData))) min(max(max(app.movieData)))], 'Parent', app.AxesStack);
            app.curSlice = s;
            app.AxesStack.Title.String = regexprep(app.curStak, '_', ' ');
            hold(app.AxesPlot, 'on');
            hTime = plot(app.AxesPlot, zeros(2,1), app.AxesPlot.YLim, 'b');
            hLeg = get(app.AxesPlot, 'Legend');
            if ~isempty(hLeg)
                hLeg.String = {hLeg.String{1}, hLeg.String{2}};
            end
            app.curTime = hTime;
            warning('on', 'all');
        end
        
        % Create a series of patch item to show the ROIs
        function showMaskDIC(app, cellRoi, varargin)
            if nargin == 3
                oldCell = varargin{1};
                cellN = app.curCell(1);
                app.patchMask(oldCell).FaceColor = [.0 .8 .8];
                app.patchMask(oldCell).EdgeColor = [.0 .8 .8];
                app.patchMask(cellN).FaceColor = [.8 .0 .8];
                app.patchMask(cellN).EdgeColor = [.8 .0 .8];
                % Show in the stack which cell we are looking at
                if app.RadioSingleTrace.Value
                    hPatch = findobj(app.AxesStack, 'Type', 'Patch');
                    if ~isempty(hPatch)
                        delete(hPatch)
                    end
                    hPatch = copyobj(app.patchMask(cellN),app.AxesStack);
                    hPatch.FaceColor = 'none';
                end
            else
                hold(app.AxesDIC, 'on')
                p = gobjects(length(cellRoi),1);
                for ij = 1:length(cellRoi)
                    p(ij) = patch(app.AxesDIC, cellRoi{ij,2}, cellRoi{ij,1}, [.0 .8 .8], 'FaceAlpha', .5, 'EdgeColor', [.0 .8 .8]);
                end
                app.patchMask = p;
                app.curCell = 1;
            end
            % Add the mask to the movie
            copyobj(app.patchMask, app.AxesStack)
            hPatch = findobj(app.AxesStack, 'Type', 'Patch');
            for p = 1:length(hPatch)
                hPatch(p).FaceColor = 'none';
            end
        end
        
        % Function to manually evaluate the drift between images
        function manualRegistrationPressed(app, event)
            whatDIC = contains(app.dicT.CellID, app.curDIC);
            dicImage2 = app.dicT.CorImage{whatDIC};
            imgFile = app.imgT.Filename{contains(app.imgT.Filename, app.curStak)};
            imgTif = imread(imgFile);
            hManReg = figure('Name', 'Manual registration');
            imshowpair(dicImage2, imgTif);
            [x, y] = ginput(2);
            app.imgT.RegEst{contains(app.imgT.Filename, app.curStak)} = [1,0,0;0,1,0;x(1)-x(2),y(1)-y(2),1];
            close(hManReg);
            ListImagesChange(app, event)
        end
        
        % Calculate intensity over the regions
        function getIntensityvalues(app)
            % First get the image where to run the analysis
            warning('off', 'all');
            togglePointer(app);
            tempMask = app.cellsMask;
            tempExp = app.dicT.ExperimentID{contains(app.dicT.CellID, app.curDIC)};
            tempRoi = app.dicT.RoiSet{contains(app.dicT.CellID, app.curDIC)};
            % Then load the each stack from this DIC
            imgFltr = find(contains(app.imgT.ExperimentID, tempExp));
            nImages = numel(imgFltr);
            hWait = waitbar(0, 'Loading image...');
            loadTime = 100;
            for i = 1:nImages
                waitbar(i/nImages, hWait, sprintf('Loading data (Remaining time ~%0.2f s)', (nImages-i+1)*loadTime));
                tic
                imgData = loadTiffFirstTime(app, imgFltr(i), hWait);
                % Adjust the orientation of the image
                %imgData = flip(imgData, 1);
                %imgData = imrotate(imgData, -90);
                nFrames = size(imgData, 3);
                % Data might need registration
                if app.options.Registration && ~contains(app.imgT{imgFltr(i), 'CellID'}, app.options.Reference)
                    if any(strcmp('RegEst', app.imgT.Properties.VariableNames))
                        if isempty(cell2mat(app.imgT{imgFltr(i),'RegEst'}))
                            tyrFile = app.imgT{find(contains(app.imgT{imgFltr,'CellID'}, app.options.Reference),1,'first'),'Filename'};
                            imgTyr = imread(cell2mat(tyrFile));
                            imgOut = imref2d(size(imgTyr));
                            imgTif = imgData(:,:,1);
                            tformEstimate = imregcorr(imgTif, imgTyr, 'translation');
                            app.imgT{imgFltr(i),'RegEst'} = {tformEstimate.T};
                        else
                            tformEstimate = affine2d(cell2mat(app.imgT{imgFltr(i),'RegEst'}));
                            imgOut = imref2d(size(imgData(:,:,1)));
                        end
                    else
                        tyrFile = app.imgT{find(contains(app.imgT{imgFltr,'CellID'}, app.options.Reference),1,'first'),'Filename'};
                        imgTyr = imread(cell2mat(tyrFile));
                        imgOut = imref2d(size(imgTyr));
                        imgTif = imgData(:,:,1);
                        tformEstimate = imregcorr(imgTif, imgTyr, 'translation');
                        app.imgT{imgFltr(i),'RegEst'} = {tformEstimate.T};
                    end
                    imgData = imwarp(imgData, tformEstimate, 'OutputView', imgOut);
                end
                % Get the Z profile of the ROIs
                nRoi = size(tempRoi,1);
                roiIntensities = zeros(nRoi, nFrames);
                for roi = 1:nRoi
                    %roiIntensities(roi,:) = mean(mean(imgData(tempMask{roi,2}(1):tempMask{roi,2}(3), tempMask{roi,1}(1):tempMask{roi,1}(2), :))); % 5 is the radius, needs to be set as option
                    roiIntensities(roi,:) = mean(mean(imgData(tempMask{roi,1}(1):tempMask{roi,1}(2), tempMask{roi,2}(1):tempMask{roi,2}(3), :))); % 5 is the radius, needs to be set as option
                end
                % Detect the minimum intensity in 10 region of the recordings to calculate the deltaF/F0
                frameDividers = [1:round(nFrames / 10):nFrames, nFrames];
                minVals = zeros(nRoi, 10);
                minIdxs = zeros(nRoi, 10);
                for idx = 1:10
                    [minVals(:,idx), minIdxs(:,idx)] = min(roiIntensities(:,frameDividers(idx):frameDividers(idx+1)), [], 2);
                end
                baseInts = mean(minVals,2);
                deltaff0Ints = (roiIntensities - repmat(baseInts, 1, nFrames)) ./ repmat(baseInts, 1, nFrames);
                app.imgT(imgFltr(i), 'RawIntensity') = {roiIntensities};
                app.imgT(imgFltr(i), 'FF0Intensity') = {deltaff0Ints};
                loadTime = toc;
            end
            detrendData(app)
            delete(hWait);
            togglePointer(app);
            warning('on', 'all');
        end
        
        % Detrend the data using a moving median or polynomial fit
        function detrendData(app)
            tempExp = app.dicT.ExperimentID{contains(app.dicT.CellID, app.curDIC)};
            %imgFltr = find(contains(app.imgT.ExperimentID, tempExp));
            imgFltr = find(~cellfun(@isempty, app.imgT.FF0Intensity));
            nImages = numel(imgFltr);
            tempDetrend = cell(nImages,1);
            for i = 1:nImages
                tempData = app.imgT{imgFltr(i), 'FF0Intensity'}{:};
                switch app.options.Detrending
                    case 'None'
                        tempDetrend{i,1} = tempData;
                    case 'Moving median'
                        fitData = movmedian(tempData, [0.2*app.options.DetrendSize 0.8*app.options.DetrendSize], 2);
                        tempDetrend{i,1} = tempData - fitData;
                    case 'Polynomial'
                        tempDetrend{i,1} = detrend(tempData', app.options.DetrendSize)';
                end
            end
            app.imgT{imgFltr, 'DetrendData'} = tempDetrend;
            app.imgT = movevars(app.imgT, 'DetrendData', 'Before', 'SpikeLocations');
        end
        
        % Create a function to load the data for the first time
        function imgData = loadTiffFirstTime(app, imgIdx, hWait)
            imgFiles = app.imgT.Filename{imgIdx};
            imgWidth = app.imgT.ImgProperties(imgIdx,1);
            imgHeight = app.imgT.ImgProperties(imgIdx,2);
            imgNumber = app.imgT.ImgProperties(imgIdx,3);
            imgData = zeros(imgWidth, imgHeight, imgNumber, 'uint16');
            %startPoint = cell2mat(app.imgT{imgIdx,'ImgByteStrip'});
            imgStackID = fopen(imgFiles, 'r');
            tstack = Tiff(imgFiles);
            imgData(:,:,1) = tstack.read();
            for n = 2:imgNumber
                %fseek(imgStackID, startPoint(n), 'bof');
                %imgData(:,:,n) = fread(imgStackID, [imgWidth, imgHeight], 'uint16=>uint16');
                nextDirectory(tstack);
                imgData(:,:,n) = tstack.read();
                %imgData(:,:,n) = imread(imgFiles,n);
            end
            fclose(imgStackID);
        end
        
        % Detect the spike in one cell
        function detectSpike(app, imgIdx)
            warning('off', 'all');
            % Gather the image data
            detrendData = cell2mat(app.imgT{imgIdx, 'DetrendData'});
            nanFrames = find(isnan(detrendData(1,:)), 1);
            if ~isempty(nanFrames)
                detrendData = detrendData(:,1:nanFrames-1);
            end
            [nTraces, nFrames] = size(detrendData);
            Fs = app.imgT.ImgProperties(imgIdx,4);
            % denoise the data using wavelet signal denoising with standard parameters
            if any(isnan(detrendData))
                detrendData(isnan(detrendData)) = 0;
            end
            try
                smoothData = wdenoise(detrendData', 'DenoisingMethod', 'BlockJS')'; % wdenoise only works on colums so inverted twice
            catch ME
                rethrow(ME)
            end
            spikeInts = cell(nTraces, 1);
            spikeLocs = cell(nTraces, 1);
            spikeWidths = cell(nTraces, 1);
            spikeRaster = nan(nTraces, nFrames);
            spikeProm = app.options.PeakMinProminance;
            spikeDist = app.options.PeakMinDistance / Fs;
            spikeMinLeng = app.options.PeakMinDuration / Fs;
            spikeMaxLeng = app.options.PeakMaxDuration / Fs;
            % for parfor settings
            switch app.options.DetectTrace
                case 'Raw'
                    gradientTrace = detrendData;
                case 'Gradient'
                    gradientTrace = gradient(smoothData);
                case 'Smooth'
                    gradientTrace = smoothData;
            end
            switch app.options.PeakMinHeight
                case 'MAD'
                    spikeThr = median(gradientTrace,2) + mad(gradientTrace,0,2) * app.options.SigmaThr;
                case 'Normalized MAD'
                    spikeThr = median(gradientTrace,2) + mad(gradientTrace,0,2) * app.options.SigmaThr * (-1 / (sqrt(2) * erfcinv(3/2)));
                case 'Rolling StDev'
                    warndlg(sprintf('Not implemented yet!\nUsed MAD instead.'), 'Failed detection')
                    spikeThr = median(gradientTrace,2) + mad(gradientTrace,0,2) * app.options.SigmaThr;
            end
            for trace = 1:nTraces
                [spikeInts{trace,1}, spikeLocs{trace,1}, spikeWidths{trace,1}] = findpeaks(gradientTrace(trace,:), Fs, 'MinPeakDistance', spikeDist, ...
                    'MinPeakHeight', spikeThr(trace), 'MinPeakProminence', spikeProm, 'MinPeakWidth', spikeMinLeng, 'MaxPeakWidth', spikeMaxLeng);
            end
            spikeRaster = spikeRaster .* repmat((1:nTraces)', 1, nFrames);
            app.imgT.SpikeLocations{imgIdx} = spikeLocs;
            app.imgT.SpikeIntensities{imgIdx} = spikeInts;
            app.imgT.SpikeWidths{imgIdx} = spikeWidths;
            app.imgT.SpikeRaster{imgIdx} = spikeRaster;
            warning('on', 'all');
            updatePlot(app);
        end
        
        % Quantification of the spike data
        function quantifySpikes(app, imgIdx)
            warning('off', 'all');
            % Create a raster plot
            Fs = app.imgT.ImgProperties(imgIdx,4);
            tempData = app.imgT.DetrendData{imgIdx};
            nFrames = size(tempData,2);
            tempLoc = cellfun(@(x) round(x * Fs), app.imgT.SpikeLocations{imgIdx}, 'UniformOutput', false);
            % Calculate the rise and decay time
            % First check if there is detection closer than the minimum distance
            diffLoc = cellfun(@diff, tempLoc, 'UniformOutput', false);
            anyLoc = cellfun(@(x) any(x <= app.options.PeakMinDistance), diffLoc);
            if any(anyLoc)
                hasSmall = find(anyLoc);
                for hS = hasSmall'
                    keepIdx = true(1,numel(app.imgT.SpikeIntensities{imgIdx}{hS}));
                    smallIdx = find(diffLoc{hS} <= app.options.PeakMinDistance);
                    [~, minIdx] = min(app.imgT.SpikeIntensities{imgIdx}{hS}(smallIdx:smallIdx+1));
                    keepIdx(minIdx+smallIdx-1) = ~keepIdx(minIdx+smallIdx-1);
                    app.imgT.SpikeLocations{imgIdx}{hS} = app.imgT.SpikeLocations{imgIdx}{hS}(keepIdx);
                    app.imgT.SpikeIntensities{imgIdx}{hS} = app.imgT.SpikeIntensities{imgIdx}{hS}(keepIdx);
                    app.imgT.SpikeWidths{imgIdx}{hS} = app.imgT.SpikeWidths{imgIdx}{hS}(keepIdx);
                end
                tempLoc = cellfun(@(x) round(x * Fs), app.imgT.SpikeLocations{imgIdx}, 'UniformOutput', false);
            end
            [spikeStart, spikeEnd, spikePeak, riseTime, decayTau, spikeLoc] = riseAndDecay(app, tempData, tempLoc, Fs);
            tempLoc = spikeLoc';
            %app.imgT.SpikeLocations{imgIdx} = cellfun(@(x) x / Fs, spikeLoc, 'UniformOutput', false);
            % Create the raster plot
            nCell = size(tempLoc, 1);
            tempRast = nan(nCell, nFrames);
            for c = 1:nCell
                sStart = spikeStart{c};
                sEnd = spikeEnd{c};
                if ~isempty(sStart)
                    for s = 1:size(sStart,2)
                        tempRast(c,sStart(2,s):min(sEnd(:,s))) = c;
                    end
                end
            end
            % Calculate the spike frequency for single cell
            totTime = (nFrames-1) / Fs;
            tempFreq = cellfun(@(x) numel(x) / totTime, tempLoc);
            % Calculate the network event with gaussian smoothing of the sum of the raster (similar to EvA)
            networkRaster = tempRast;
            networkRaster(~isnan(networkRaster)) = 1;
            networkRaster(isnan(networkRaster)) = 0;
            networkRaster = sum(networkRaster);
            smoothWindow = gausswin(10);
            smoothWindow = smoothWindow / sum(smoothWindow);
            networkRaster = filter(smoothWindow, 1, networkRaster);
            [networkPeaks, networkLocs, networkFWHM] = findpeaks(networkRaster, Fs, 'WidthReference', 'halfheight');
            % Calculate the synchronous frequency (cells firing +- 1 frames) ### NEEDS UPDATING ###
            allLocs = sort(cell2mat(tempLoc'));
            uniLocs = unique(allLocs);
            synLocs = [];
            synPC = [];
            synLoc = 1;
            for loc = 1:numel(uniLocs)
                locRange = [];
                if loc == 1
                    locRange = [uniLocs(loc), uniLocs(loc)+2];
                elseif loc == numel(uniLocs)
                    locRange = [uniLocs(loc)-2, uniLocs(loc)];
                else
                    if uniLocs(loc) - uniLocs(loc-1) > 1
                        locRange = [uniLocs(loc), uniLocs(loc)+2];
                    end
                end
                if ~isempty(locRange)
                    synLocs(synLoc) = uniLocs(loc)/Fs;
                    synPC(synLoc) = sum(allLocs>=locRange(1) & allLocs<=locRange(2))/nCell*100;
                    synLoc = synLoc + 1;
                end
            end
            % Calculate the network frequency as the number of spikes that have > 80% of neurons firing
            netThreshold = sum(any(~isnan(tempRast), 2)) * 0.8;
            networkFreq = sum(synPC >= netThreshold) / totTime;
            networkFreqGauss = sum(networkPeaks >= netThreshold) / totTime; % Use the maximum number of cells
            % Calculate the interspike intervals, and time-frequency
            interSpikeInterval = cellfun(@(x) diff(x) / Fs, tempLoc, 'UniformOutput', false);
            qcd = @(x) (quantile(x, .75) - quantile(x, .25)) / (quantile(x, .75) + quantile(x, .25));
            ISI_CoV = cellfun(qcd, interSpikeInterval);
            nSample = 5;
            timeISI = 0:nFrames/nSample:nFrames;
            tempMeanISI = nan(size(tempLoc,1),nSample);
            tempMeanFreq = nan(size(tempLoc,1),nSample);
            for s = 1:size(tempLoc,1)
                if ~isempty(tempLoc{s})
                    tempISI = [tempLoc{s}(1) / Fs, interSpikeInterval{s}];
                    for t = 1:nSample
                        tempMeanISI(s,t) = mean(tempISI(tempLoc{s} > timeISI(t) & tempLoc{s} < timeISI(t+1)));
                    end
                end
            end
            % Calculate and keep track of subthreshold spikes (as for now defined as promininence < 0.2 a.u.)
            subFltr = cellfun(@(x) x<=0.2, spikePeak, 'UniformOutput', false);
            subRatio = cellfun(@(x) sum(x)/numel(x), subFltr);
            % Calculate if there are bursts
            
            % Save the data to the table
            app.imgT.SpikeRaster{imgIdx} = tempRast;
            app.imgT.Participation(imgIdx) = sum(any(~isnan(tempRast), 2)) / nCell * 100;
            app.imgT.CellFrequency{imgIdx} = tempFreq;
            app.imgT.SubthresholdFltr{imgIdx} = subFltr';
            app.imgT.SubthresholdRatio{imgIdx} = subRatio';
%             app.imgT.SynchronousLocations{imgIdx} = synLocs;
%             app.imgT.SynchronousPercentages{imgIdx} = synPC;
            app.imgT.NetworkFrequency(imgIdx) = networkFreqGauss;
            app.imgT.InterSpikeInterval{imgIdx} = interSpikeInterval;
            app.imgT.InterSpikeIntervalCoV{imgIdx} = ISI_CoV;
%             app.imgT.TimeInterSpikeInterval{imgIdx} = nanmean(tempMeanISI);
            app.imgT.NetworkPeaks{imgIdx} = networkPeaks;
            app.imgT.NetworkLocs{imgIdx} = networkLocs;
%             app.imgT.NetworkFWHMGauss{imgIdx} = networkFWHM;
            app.imgT.NetworkRaster{imgIdx} = networkRaster;
            app.imgT.SpikeStartAndEnd{imgIdx} = cellfun(@(x, y) [x;y], spikeStart, spikeEnd, 'UniformOutput', false);
            app.imgT.SpikeWidths{imgIdx} = cellfun(@(x,y) (y(2,:) - x(2,:)) / Fs, spikeStart, spikeEnd, 'UniformOutput', false)';
            app.imgT.SpikeIntensities{imgIdx} = spikePeak';
            app.imgT.CellRiseTime{imgIdx} = riseTime';
            app.imgT.CellDecayCon{imgIdx} = cellfun(@(x) x(1,:), decayTau, 'UniformOutput', false)';
            app.imgT.CellDecayTau{imgIdx} = cellfun(@(x) x(2,:), decayTau, 'UniformOutput', false)';
            warning('on', 'all');
        end
        
        % Define rise and decay time
        function [spikeStart, spikeEnd, spikePeak, spikeRise, spikeDecay, spikeLoc] = riseAndDecay(app, traceData, spikeData, Fs)
            % Loop through all the identified peaks to find the left and right bounds
            nTraces = size(traceData,1);
            minD = app.options.PeakMinDuration;
            spikeStart = cell(1, nTraces);
            spikeLoc = cell(1, nTraces);
            spikeEnd = cell(1, nTraces);
            spikePeak = cell(1, nTraces);
            spikeRise = cell(1, nTraces);
            spikeDecay = cell(1, nTraces);
            for t = 1:nTraces
                smoothData = wdenoise(traceData(t,:), 'DenoisingMethod', 'BlockJS');
                % Calculate all the local minima
                allValleys = islocalmin(smoothData,2);
                spikeLocs = spikeData{t} + 1;
                nSpikes = numel(spikeLocs);
                indexLB = nan(2, nSpikes);
                indexRB = nan(2, nSpikes);
                promInt = nan(1, nSpikes);
                decayTauA = nan(1, nSpikes);
                decayTauB = nan(1, nSpikes);
                for s = 1:nSpikes
                    % Might be useful to refine to spike location
                    [newMax, newSpikeLoc] = max(traceData(t,max(1,spikeLocs(s)-minD):min(length(traceData),spikeLocs(s)+10)));
                    if newMax > traceData(t,spikeLocs(s))
                        if s < nSpikes && (newSpikeLoc + spikeLocs(s) - minD - 1) < (spikeLocs(s+1)-2)
                            spikeLocs(s) = newSpikeLoc + spikeLocs(s) - minD - 1;
                        end
                    end
                    % Left bound
                    if s == 1
                        tempStart = 1;
                    else
                        tempStart = spikeLocs(s-1);
                    end
                    tempIdx = find(allValleys(tempStart:spikeLocs(s)), 1, 'last')  + tempStart;
                    if isempty(tempIdx)
                        tempIdx = tempStart;
                    elseif(spikeLocs(s) - tempIdx <= minD)
                        tempIdx = find(allValleys(tempStart:tempIdx-minD), 1, 'last') + tempStart;
                        if isempty(tempIdx)
                            tempIdx = min(tempStart + app.options.PeakMinDistance, find(allValleys(tempStart:spikeLocs(s)), 1, 'last')  + tempStart);
                        end
                    end
                    indexLB(1,s) = tempIdx;
                    tempHalfMax = mean([smoothData(tempIdx), smoothData(spikeLocs(s))]);
                    % Right bound
                    if s == nSpikes
                        tempEnd = length(smoothData);
                    else
                        tempEnd = spikeLocs(s+1);
                    end
                    tempIdx = find(allValleys(spikeLocs(s):tempEnd), 1, 'first') + spikeLocs(s);
                    if isempty(tempIdx)
                        tempIdx = tempEnd;
                    else
                        while (tempIdx - spikeLocs(s) <= minD * 2) || (smoothData(tempIdx) > tempHalfMax)
                            tempIdx = find(allValleys((tempIdx+minD*2):tempEnd), 1, 'first') + tempIdx+(minD*2);
                            if isempty(tempIdx)
                                tempIdx = max(tempEnd - app.options.PeakMinDistance, find(allValleys(spikeLocs(s):tempEnd), 1, 'first') + spikeLocs(s));
                                break
                            end
                        end
                    end
                    indexRB(1,s) = tempIdx;
                    % Calculate the FWHM (based on the prominence)
                    baseInt = min(traceData(t, indexLB(1,s):indexRB(1,s)));
                    promInt(s) = traceData(t, spikeLocs(s)) - baseInt;
                    halfMax = baseInt + promInt(s)/2;
                    tempIdx = find(smoothData(indexLB(1,s):indexRB(1,s)) > halfMax, 1, 'first') - 2 + indexLB(1,s);
                    if isempty(tempIdx) || tempIdx == 0
                        tempIdx = 1;
                    end
                    indexLB(2,s) = tempIdx;
                    tempIdx = find(smoothData(indexLB(1,s):indexRB(1,s)) > halfMax, 1, 'last') + 2 + indexLB(1,s);
                    if isempty(tempIdx)
                        tempIdx = length(smoothData);
                    end
                    indexRB(2,s) = tempIdx;
                end
                % Work with parallel for loop to calculate the decay fit
                tempTrace = traceData(t,:);
                parfor s = 1:nSpikes
                    % Calculate the decay constant
                    decayTime = (spikeLocs(s):indexRB(1,s)) / Fs;
                    decayTrace = tempTrace(spikeLocs(s):indexRB(1,s));
                    expFit = fittype( 'exp1' );
                    fitOpts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    fitOpts.Display = 'Off';
                    fitOpts.Normalize = 'On';
                    fitOpts.StartPoint = [0 0];
                    if numel(decayTrace) > 2
                        decayFit = fit(decayTime', decayTrace', expFit, fitOpts);
                        % To get the fitted trace use: fitTrace = decayFit.a * exp(decayFit.b .* ((decayTime-mean(decayTime))/std(decayTime)));
                        decayTauA(1,s) = decayFit.a;
                        decayTauB(1,s) = -1*(decayFit.b)^-1;
                    end
                end
                spikeStart{t} = indexLB;
                spikeLoc{t} = spikeLocs;
                spikeEnd{t} = indexRB;
                spikePeak{t} = promInt;
                spikeRise{t} = (spikeLocs - indexLB(1,:)) / Fs;
                spikeDecay{t} = [decayTauA; decayTauB];
            end
        end
            
        % Update the plot
        function updatePlot(app)
            % First get the data
            imgID = contains(app.imgT.Filename, app.curStak);
            tempData = app.imgT.DetrendData{imgID};
            if isempty(tempData)
                return
            end
            switch app.options.DetectTrace
                case 'Gradient'
                    tempData = gradient(tempData);
                case 'Smooth'
                    tempData = wdenoise(tempData', 'DenoisingMethod', 'BlockJS')';
            end
            Fs = app.imgT.ImgProperties(imgID,4);
            time = (0:size(tempData, 2)-1) / Fs;
            if app.RadioAllMean.Value == true
                hold(app.AxesPlot, 'off')
                hLeg(1) = plot(app.AxesPlot, time, tempData(1,:), 'Color', [.7 .7 .7]);
                hold(app.AxesPlot, 'on')
                plot(app.AxesPlot, time, tempData(2:end,:), 'Color', [.7 .7 .7]);
                hLeg(2) = plot(app.AxesPlot, time, nanmean(tempData), 'Color', 'r');
                legend(app.AxesPlot, hLeg, {'All' 'Mean'}, 'Location', 'best', 'box', 'off')
                app.AxesPlot.XLim = [time(1) time(end)];
                app.AxesPlot.Title.String = 'Ca^{2+} Traces';
                app.AxesPlot.YLabel.String = '\DeltaF/F_0';
                app.AxesPlot.XLabel.String = 'Time (s)';
                app.AxesPlot.TickDir = 'out';
                app.AxesPlot.Box = 'off';
            end
            if app.RadioSingleTrace.Value == true
                % plot one trace and the identified spikes
                cellN = app.curCell(1);
                hold(app.AxesPlot, 'off')
                plot(app.AxesPlot, time, tempData(cellN,:), 'Color', 'k');
                hold(app.AxesPlot, 'on')
                if ~isempty(app.imgT.SpikeLocations{imgID})
                    spikeLoc = app.imgT.SpikeLocations{imgID}{cellN};
                    spikeInt = tempData(cellN, round(spikeLoc*Fs)+1);
                    % spikeInt = app.imgT.SpikeIntensities{imgID}{cellN};
                    plot(app.AxesPlot, spikeLoc, spikeInt, 'or')
                end
                % Plot the rise and decay time if calculated
                if any(contains(app.imgT.Properties.VariableNames, 'SpikeStartAndEnd')) && ~isempty(app.imgT.SpikeStartAndEnd{imgID}{cellN})
                    spikeStart = app.imgT.SpikeStartAndEnd{imgID}{cellN}(1,:);
                    spikeEnd = app.imgT.SpikeStartAndEnd{imgID}{cellN}(3,:);
                    nSpike = numel(spikeLoc);
                    for s = 1:nSpike
                        spikeRiseTime = spikeStart(s):round(spikeLoc(s)*Fs)+1;
                        spikeDecayTime = round(spikeLoc(s)*Fs)+1:spikeEnd(s);
                        plot(app.AxesPlot, time(spikeRiseTime), tempData(cellN,spikeRiseTime), 'm');
                        plot(app.AxesPlot, time(spikeDecayTime), tempData(cellN,spikeDecayTime), 'g');
                    end
                end
                legend(app.AxesPlot, 'hide')
                % Add the level of the threshold
                switch app.options.PeakMinHeight
                    case 'MAD'
                        spikeThr = median(tempData(cellN,:)) + mad(tempData(cellN,:)) * app.options.SigmaThr;
                    case 'Normalized MAD'
                        spikeThr = median(tempData(cellN,:)) + mad(tempData(cellN,:)) * app.options.SigmaThr * (-1 / (sqrt(2) * erfcinv(3/2)));
                    case 'Rolling StDev'
                        warndlg(sprintf('Not implemented yet!\nUsed MAD instead.'), 'Failed detection')
                        spikeThr = median(tempData(cellN,:)) + mad(tempData(cellN,:)) * app.options.SigmaThr;
                end
                plot(app.AxesPlot, app.AxesPlot.XLim, [spikeThr spikeThr], ':b')
                % show in the DIC which one whe are looking at
                if numel(app.curCell) == 2
                    oldCell = app.curCell(2);
                else
                    oldCell = 1;
                end
                showMaskDIC(app, app.cellsMask, oldCell)
            end
            app.AxesPlot.XLim = [time(1) time(end)];
            if app.ToggleFixAxis.Value
                app.AxesPlot.YLim = app.YMinMax;
            end
            app.AxesPlot.Title.String = 'Ca^{2+} Traces';
            app.AxesPlot.YLabel.String = '\DeltaF/F_0';
            app.AxesPlot.XLabel.String = 'Time (s)';
            app.AxesPlot.TickDir = 'out';
            app.AxesPlot.Box = 'off';
            % Check if the stack is visible, and if so add the data line
            if app.ToggleViewStack.Value
                hTime = plot(app.AxesPlot, ones(2,1)*app.SliderImageStack.Value, app.AxesPlot.YLim, 'b');
                app.curTime = hTime;
                hLeg = get(app.AxesPlot, 'Legend');
                hLeg.String = {hLeg.String{1}, hLeg.String{2}};
            end
            app.AxesPlot.Toolbar.Visible = 'on';
        end
        
        % Brightness and contrast (similar to ImageJ)
        function adjustBrightness(app, img1, img2)
            % First get the histogram
            [dicCount, dicLoc] = imhist(img1);
            [stackCount, stackLoc] = imhist(img2);
            bcFig = figure('Units', 'pixels', 'Visible', 'on',...
                'Position', [500 500 403 280], 'MenuBar', 'none');
            uicontrol(bcFig, 'Style', 'text', 'String', 'DIC',...
                'FontSize', 15, 'Position', [60 260 100 20]);
                hHistDIC = axes('Units', 'pixels', 'Position', [5 150 190 100]);
                uicontrol(bcFig, 'Style', 'text', 'String', 'Minimum',...
                    'Position', [5 110 190 20]);
                hMinDIC = uicontrol(bcFig, 'Style', 'slider', 'Units', 'pixels',...
                    'Position', [5 95 190 15], 'Value', dicLoc(1), 'SliderStep', [1/numel(dicLoc) 1/numel(dicLoc)],...
                    'Min', dicLoc(1), 'Max', dicLoc(end));
                uicontrol(bcFig, 'Style', 'text', 'String', 'Maximum',...
                    'Position', [5 70 190 20]);
                hMaxDIC = uicontrol(bcFig, 'Style', 'slider', 'Units', 'pixels',...
                    'Position', [5 55 190 15], 'Value', dicLoc(end), 'SliderStep', [1/numel(dicLoc) 1/numel(dicLoc)],...
                    'Min', dicLoc(1), 'Max', dicLoc(end));
            uicontrol(bcFig, 'Style', 'text', 'String', 'Movie',...
                'FontSize', 15, 'Position', [260 260 100 20]);
                hHistStack = axes('Units', 'pixels', 'Position', [210 150 190 100]);
                uicontrol(bcFig, 'Style', 'text', 'String', 'Minimum',...
                    'Position', [210 110 190 20]);
                hMinStack = uicontrol(bcFig, 'Style', 'slider', 'Units', 'pixels',...
                    'Position', [210 95 190 15], 'Value', stackLoc(1), 'SliderStep', [1/numel(stackLoc) 1/numel(stackLoc)],...
                    'Min', stackLoc(1), 'Max', stackLoc(end));
                uicontrol(bcFig, 'Style', 'text', 'String', 'Maximum',...
                    'Position', [210 70 190 20]);
                hMaxStack = uicontrol(bcFig, 'Style', 'slider', 'Units', 'pixels',...
                    'Position', [210 55 190 15], 'Value', stackLoc(end), 'SliderStep', [1/numel(stackLoc) 1/numel(stackLoc)],...
                    'Min', stackLoc(1), 'Max', stackLoc(end));
            % DIC plot
            stairs(hHistDIC,dicLoc, dicCount, 'k', 'LineWidth', 2);
            xlim(hHistDIC, [dicLoc(1) dicLoc(end)])
            hold(hHistDIC, 'on')
            hHistDIC.YLim = [min(dicCount) max(dicCount)];
            hHistDIC.XTick = [];
            hHistDIC.YTick = [];
            hLineDIC = plot(hHistDIC, [dicLoc(1) dicLoc(end)], [min(dicCount) max(dicCount)], 'b');
            hMinDIC.Callback = {@hMovedDIC};
            hMaxDIC.Callback = {@hMovedDIC};
            addlistener(hMinDIC,'ContinuousValueChange',@(src,event)hMovedDIC(src,event));
            addlistener(hMaxDIC,'ContinuousValueChange',@(src,event)hMovedDIC(src,event));
            % Movie plot
            stairs(hHistStack, stackLoc, stackCount, 'k', 'LineWidth', 2);
            xlim(hHistStack, [stackLoc(1) stackLoc(end)])
            hHistStack.YLim = [min(stackCount) max(stackCount)];
            hold(hHistStack, 'on')
            hHistStack.XTick = [];
            hHistStack.YTick = [];
            hLineStack = plot(hHistStack, [stackLoc(1) stackLoc(end)], [min(stackCount) max(stackCount)], 'b');
            hMinStack.Callback = {@hMovedStack};
            hMaxStack.Callback = {@hMovedStack};
            addlistener(hMinStack,'ContinuousValueChange',@(src,event)hMovedStack(src,event));
            addlistener(hMaxStack,'ContinuousValueChange',@(src,event)hMovedStack(src,event));
            % Functions
            function hMovedDIC(src, event)
                app.hDIC.CData = imadjust(img1, [hMinDIC.Value hMaxDIC.Value]/2^16);
                hLineDIC.XData = [hMinDIC.Value hMaxDIC.Value];
            end
            function hMovedStack(src, event)
                app.hStack.CData = imadjust(img2, [hMinStack.Value hMaxStack.Value]/2^16);
                hLineStack.XData = [hMinStack.Value hMaxStack.Value];
            end
        end
        
        % Series of plotting functions for data visualization
        function rasterPlot(app, varX, varY, varN)
            figure()
            rAx = axes();
            hold on
            area(rAx, varX, varN, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', .3)
            plot(rAx, varX, varY, 'k', 'LineWidth', 2)
            plot(rAx, [varX(1) varX(end)], ones(2,1) * sum(any(~isnan(varY), 2)) * 0.8, 'r')
            rAx.Box = 'off';
            rAx.YLim = [0, size(varY,1)];
            rAx.YTick = 1:size(varY,1);
            rAx.TickDir = 'out';
            rAx.YLabel.String = 'Cell number';
            rAx.XLabel.String = 'Time (s)';
            rAx.Title.String = regexprep(app.curStak, '_', ' ');
        end
   
        % Plot for ISI over tiome
        function PlotTimeFrequencySelected(app, event)
            % get the data
            Fs = app.imgT.ImgProperties(1,4);
            nFrames = size(app.imgT.DetrendData{1},2);
            nSample = 5;
            timeISI = 0:nFrames/nSample:nFrames;
            semX = [timeISI(1:end-1);timeISI(2:end);timeISI(2:end);timeISI(1:end-1)] / Fs;
            conditions = categorical(app.imgT.Condition);
            uniCond = categories(conditions);
            nCond = numel(uniCond);
            allData = cell2mat(app.imgT.TimeInterSpikeInterval);
%             allData(isnan(allData)) = 60;
            figure; hAx1 = subplot(2,1,1); hold on
            hAx2 = subplot(2,1,2); hold on
            cmap = lines;
            for c = 1:nCond
                % Plot the ISI
                tempData = allData(conditions == uniCond(c), :);
                tempMean = nanmean(tempData);
                tempSEM = nanstd(tempData) / sqrt(sum(conditions == uniCond(c)));
                semY = [(tempMean-tempSEM); (tempMean-tempSEM); (tempMean+tempSEM); (tempMean+tempSEM)];
                patch(hAx1, semX, semY, cmap(c,:), 'EdgeColor', 'none', 'FaceAlpha', .3);
                hLeg1(c) = stairs(hAx1, timeISI / Fs, [tempMean, tempMean(end)], 'Color', cmap(c,:));
                % From the ISI calculate an average frequency and plot it
                tempMean = nanmean(1 ./ tempData);
                tempSEM = nanstd(1 ./ tempData) / sqrt(sum(conditions == uniCond(c)));
                semY = [(tempMean-tempSEM); (tempMean-tempSEM); (tempMean+tempSEM); (tempMean+tempSEM)];
                patch(hAx2, semX, semY, cmap(c,:), 'EdgeColor', 'none', 'FaceAlpha', .3);
                stairs(hAx2, timeISI / Fs, [tempMean, tempMean(end)], 'Color', cmap(c,:));
            end
            set(hAx1, 'TickDir', 'out')
            set(hAx2, 'TickDir', 'out')
            ylabel(hAx1, 'Inter spike interval (s)')
            ylabel(hAx2, 'Frequency (Hz)')
            xlabel(hAx2, 'Time (s)')
            legend(hAx1, hLeg1, uniCond, 'Box', 'off', 'Location', 'best')
        end
        
        % Plot a general overview of the data and the coefficient of variation
        function PlotOverviewSelected(app, event)
            % Get the conditions
            myCondition = categorical(app.imgT.Condition);
            uniCond = categories(myCondition);
            nCond = numel(uniCond);
            % Check if there are multiple recordings per FOV
            app.imgT.RecID = categorical(app.imgT.RecID);
            recID = categories(app.imgT.RecID);
            nRec = numel(recID);
            % Get the replicas
            batchList = app.imgT.Week;
            weeks = unique(batchList);
            nWeeks = numel(weeks);
            % Ask wich condition is the control
            if ~isfield(app.options, 'ConditionOrder')
                uniCond = unique(app.imgT.Condition);
                promptCond = 'Enter space-separated number for conditions:/n';
                for c=1:numel(uniCond)
                    promptCond = sprintf('%s, %s', promptCond, uniCond{c});
                end
                condOrd = inputdlg(promptCond, 'Order conditions');
                app.options.ConditionOrder = condOrd;
            end
            [contIdx, ~] = listdlg('ListString', uniCond, 'PromptString', 'Select control condition', 'SelectionMode', 'single');
            contID = uniCond{contIdx};
            myCondition = reordercats(myCondition, str2num(app.options.ConditionOrder{1}));
            
            %%% CONDITIONAL FINGERPRINT%%%
            % Select the variables and create a list of labels
            varID = [17 21 31:34 39:42 47 49 50 52 53];
            varNames = {'Participation','Network Frequency','Intensity','Duration','Rise Time','Decay \tau','Sub Intensity','Sub Duration','Sub Rise Time','Sub Decay \tau','Inter spike interval','ISI CoV','Cell Frequency','# of Cell','Peak Participation'};
            nFeature = numel(varID);
            % Define the dimensions and create the space holders
            matValue = [];
            groupNames = cell(1,(nRec*(nCond-1)));
            gr=1;
            for r = 1:nRec
                recT = app.imgT(app.imgT.RecID == recID(r), :);
                recT.Condition = myCondition(app.imgT.RecID == recID(r));
                recT.Condition = removecats(recT.Condition);
                conds = unique(recT.Condition);
                nCond = numel(conds);
                tempValue = zeros(nCond-1, nFeature);
                % Calculate the matrix of differences
                for v = 1:nFeature
                    tempData = recT{:,varID(v)};
                    % Loop through the weeks to normalize the data per culture batch
                    for w = 1:nWeeks
                        weekFltr = recT.Week == weeks(w);
                        controlFltr = recT.Condition == contID;
                        tempMean = mean(tempData(controlFltr & weekFltr), 'omitnan');
                        if tempMean > 0 % required in the case that there is no network activity (as in 1mM)
                            tempData(weekFltr) = tempData(weekFltr) / tempMean;
                        end
                    end
                    % Store the data in a matrix
                    for c = 2:nCond
                        tempValue(c-1, v) = log2(mean(tempData(recT.Condition == conds(c) & tempData > 0), 'omitnan'));
                        if v == 1
                            groupNames{gr} = sprintf('%s %s', conds(c), recID{r});
                            gr = gr+1;
                        end
                    end
                end
                matValue = [matValue; tempValue];
            end
            % Plot the conditional fingerprints
            % Create divergin colormap, from [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
            m1 = 10000*0.5;
            r = (0:m1-1)'/max(m1-1,1);
            g = r;
            r = [r; ones(m1,1)];
            g = [g; flipud(g)];
            b = flipud(r);
            cmap = flipud([r g b]);
            figure('Name', 'Conditional fingeprint')
            colormap(cmap)
            maxVal = max(abs([min(min(matValue)), max(max(matValue))]));
            imagesc([-repmat(maxVal, (nCond-1)*nRec, 1) matValue repmat(maxVal, (nCond-1)*nRec, 1)]);
            xlim([1.5 numel(varID)+1.5])
            colorbar
            set(gca, 'XTick', 2:numel(varID)+1)
            set(gca, 'XTickLabel', varNames)
            set(gca, 'XTickLabelRotation', 45)
            set(gca, 'YTick', 1:(nCond-1)*nRec)
            set(gca, 'YTickLabel', groupNames)
            box off
            set(gca, 'TickDir', 'out')
            hold on
            for r = 0.5:1:((nCond-1)*nRec)+0.5
                plot([0 nFeature+2], [r r], 'k')
            end
            for c = 1.5:1:nFeature+1.5
                plot([c c], [0 nCond*nRec], 'k')
            end
            
            %%%BOXPLOTS OF FINGERPRINT%%%
            % Calculate the 95% confidence interval (of the mean) and overlay to the graph
            CIFcn = @(x,p) nanstd(x)/sqrt(sum(~isnan(x))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x))-1) + nanmean(x);
            for r = 1:nRec
                recT = app.imgT(app.imgT.RecID == recID(r), :);
                recT.Condition = myCondition(app.imgT.RecID == recID(r));
                recT.Condition = removecats(recT.Condition);
                uniCond = categories(recT.Condition);
                nCond = numel(uniCond);
                figure('Name', sprintf('Overview of recording: %s', recID{r}))
                cmap = lines;
                for p = 1:nFeature
                    % Create the main plot
                    subplot(3,5,p); hold on
                    yData = recT{:,varID(p)};
%                     QCD = nan(nCond, nWeeks);
                    % CI
                    CI_Cont = CIFcn(recT{recT.Condition==contID, varID(p)}, 95);
                    patch([0 nCond+1 nCond+1 0], [CI_Cont(1) CI_Cont(1) CI_Cont(2) CI_Cont(2)], cmap(1,:), 'EdgeColor', 'none', 'FaceAlpha',.1)
                    for c = 1:nCond
                        condFltr = recT.Condition == uniCond(c);
                        tempY = sort(yData(condFltr));
                        quantY = quantile(tempY, [0.25 0.5 0.75]);
                        minW = quantY(1) - 1.5*(quantY(3)-quantY(1));
                        lowW = find(tempY>=minW,1,'first');
                        minW = tempY(lowW);
                        maxW = quantY(3) + 1.5*(quantY(3)-quantY(1));
                        highW = find(tempY<=maxW,1,'last');
                        maxW = tempY(highW);
                        % Boxplot
                        patch([c-.25 c+.25 c+.25 c-.25], [quantY(1) quantY(1) quantY(3) quantY(3)], cmap(c,:), 'FaceAlpha', .3, 'EdgeColor', cmap(c,:));
                        plot([c-.25 c+.25], [quantY(2) quantY(2)], 'color', cmap(c,:), 'LineWidth', 2);
                        plot([c c], [minW quantY(1)], 'color', cmap(c,:));
                        plot([c c], [quantY(3) maxW], 'color', cmap(c,:));
                        % Add the data points
                        x = linspace(c - 0.15, c + 0.15, nWeeks);
                        for w = 1:nWeeks
                            weekFltr = recT.Week == weeks(w);
                            if sum(weekFltr & condFltr) > 0
                                plot(x(w),yData(weekFltr & condFltr), 'o', 'MarkerEdgeColor', cmap(c,:), 'MarkerSize',4,'MarkerFaceColor','w')
                                % calculate the quartile coefficient of dispersion
%                                 QCD(c,w) = (quantile(yData(weekFltr & condFltr), 0.75) - quantile(yData(weekFltr & condFltr), 0.25)) / (quantile(yData(weekFltr & condFltr), 0.75) + quantile(yData(weekFltr & condFltr), 0.25)) * 100;
                            end
                        end
                    end
                    % Add the label
                    set(gca, 'TickDir', 'out');
                    xlim([.5 nCond+.5])
                    set(gca, 'XTick', 1:nCond);
                    set(gca, 'XTickLabel', uniCond);
                    set(gca, 'XTickLabelRotation', 45)
                    ylabel(varNames(p))
%                     % Add the coefficent of variation (or better quartile coefficient of dispersion) calculated per week
%                     subplot(2, nPlot, p+nPlot); hold on
%                     for c = 1:nCond
%                         tempY = sort(QCD(c,:));
%                         quantY = quantile(tempY, [0.25 0.5 0.75]);
%                         minW = quantY(1) - 1.5*(quantY(3)-quantY(1));
%                         lowW = find(tempY>=minW,1,'first');
%                         minW = tempY(lowW);
%                         maxW = quantY(3) + 1.5*(quantY(3)-quantY(1));
%                         highW = find(tempY<=maxW,1,'last');
%                         maxW = tempY(highW);
%                         % Boxplot
%                         patch([c-.25 c+.25 c+.25 c-.25], [quantY(1) quantY(1) quantY(3) quantY(3)], cmap(c,:), 'FaceAlpha', .3, 'EdgeColor', cmap(c,:));
%                         plot([c-.25 c+.25], [quantY(2) quantY(2)], 'color', cmap(c,:), 'LineWidth', 2);
%                         plot([c c], [minW quantY(1)], 'color', cmap(c,:));
%                         plot([c c], [quantY(3) maxW], 'color', cmap(c,:));
%                         % Add the data points
%                         x = linspace(c - 0.15, c + 0.15, nWeeks);
%                         for w = 1:nWeeks
%                             weekFltr = batchList == weeks(w);
%                             if sum(weekFltr & condFltr) > 0
%                                 plot(x(w),QCD(c,w), 'o', 'MarkerEdgeColor', cmap(c,:), 'MarkerSize',4,'MarkerFaceColor','w')
%                             end
%                         end
%                     end
%                     % Add the label
%                     set(gca, 'TickDir', 'out');
%                     xlim([.5 nCond+.5])
%                     set(gca, 'XTick', 1:nCond);
%                     set(gca, 'XTickLabel', uniCond);
%                     ylabel(yLabels(2,p))
                end
            end
        end 
    end
    
    % Callbacks methods
    methods (Access = private)
        function FileMenuOpenSelected(app, event)
            % Locate the folder with the images
            hOpt = OptionMenuSelected(app, event.Source.Text);
            waitfor(hOpt);
            imgPath = uigetdir(app.options.LastPath, 'Select Image folder');
            togglePointer(app)
            if imgPath ~= 0
                % Save the path to the settings
                app.options.LastPath = imgPath;
                % Load the data starting from the DIC/BF/Still image
                hWait = waitbar(0, 'Loading images data');
                app.DicMenu.Enable = 'on';
                imgFiles = dir(fullfile(imgPath, '*.tif'));
                % Populate the DIC table
                dicFltr = find(contains({imgFiles.name}, app.options.StillName))';
                nameParts = regexp({imgFiles.name}, '_', 'split')';
                if isempty(dicFltr)
                    % It might be that there is no a DIC image. Ask the user for the name of the image to detect cells.
                    % For this consider a format or the type YYMMDD_Condition_cs_rec
                    conditionIDs = unique(cellfun(@(x) x(2), nameParts));
                    [condIdx, OKed] = listdlg('PromptString', 'Select condition of still image', 'ListString',conditionIDs);
                    if OKed
                        app.options.StillName = conditionIDs{condIdx};
                        dicFltr = find(contains({imgFiles.name}, app.options.StillName));
                    else
                        return
                    end
                    
                end
                tempT = table;
                tempT.Filename = fullfile({imgFiles(dicFltr).folder}', {imgFiles(dicFltr).name}');
                tempT.CellID = cellfun(@(x) x(1:end-4), {imgFiles(dicFltr).name}', 'UniformOutput', false);
                if strcmp(app.options.StillName, 'DIC.')
                    expIDs = cellfun(@(x) sprintf('%s_%s_%s', x{1}, x{2}, x{3}), nameParts, 'UniformOutput', false);
                    tempT.ExperimentID = expIDs(dicFltr);
                else
%                     expIDs = cellfun(@(x) sprintf('%s_%s_%s', x{1}, x{3}, x{4}), nameParts, 'UniformOutput', false);
                    expIDs = cellfun(@(x) sprintf('%s_%s', x{1}, x{3}), nameParts, 'UniformOutput', false);
                    tempT.ExperimentID = expIDs(dicFltr);
                end
                tempDicImages = cell(size(tempT, 1), 2);
                for i = 1:size(tempT,1)
                    waitbar(i/numel(dicFltr), hWait, sprintf('Loading DIC data %0.2f%%', i/numel(dicFltr)*100));
                    dicFile = tempT.Filename{i};
                    dicImage = imread(dicFile);
                    % Correct DIC image for bad illumination
                    dicImage2 = imcomplement(dicImage);
                    filterSE = strel('square', 10);
                    dicBackground = imerode(dicImage2, filterSE);
                    dicBackground = imdilate(dicBackground, filterSE);
                    dicBackground = imgaussfilt(dicBackground, 10);
                    dicBackground = dicBackground - min(min(dicBackground));
                    dicImage2 = imcomplement(dicImage2 - dicBackground);
                    tempDicImages{i,1} = dicImage;
                    tempDicImages{i,2} = dicImage2;
                end
                tempT.RawImage = tempDicImages(:,1);
                tempT.CorImage = tempDicImages(:,2);
                tempT.RoiSet = repmat({[]}, size(tempT,1), 1);
                % Store the dicT in the app
                app.dicT = tempT;
                
                % Populate the imgT
                imgFltr = find(~contains({imgFiles.name}, app.options.StillName))';
                app.stackDs = cell(numel(imgFltr),1);
                tempT = cell(numel(imgFltr)+1, 14);
                tempT(1,:) = {'Filename', 'CellID', 'Week', 'CoverslipID', 'RecID', 'Condition', 'ExperimentID', 'ImgProperties', 'ImgByteStrip', 'RawIntensity', 'FF0Intensity', 'SpikeLocations', 'SpikeIntensities', 'SpikeWidths'};
                tempT(2:end,1) = fullfile({imgFiles(imgFltr).folder}, {imgFiles(imgFltr).name});
                tempT(2:end,2) = cellfun(@(x) x(1:end-4), {imgFiles(imgFltr).name}, 'UniformOutput', false);
                imgIDs = nameParts(imgFltr);
                for i = 1:numel(imgFltr)
                    waitbar(i/numel(imgFltr), hWait, sprintf('Loading movie data %0.2f%%', i/numel(imgFltr)*100));
                    tempT{i+1,3} = weeknum(datetime(imgIDs{i}{1}, 'InputFormat', 'yyMMdd'));
                    if strcmp(app.options.StillName, 'DIC.')
                        tempT{i+1,4} = imgIDs{i}{2};
                        tempT{i+1,5} = imgIDs{i}{3};
                        tempT{i+1,6} = imgIDs{i}{4};
                        tempT{i+1,7} = [imgIDs{i}{1} '_' imgIDs{i}{2} '_' imgIDs{i}{3}]; % use to link the DIC to the movies
                    else
                        tempT{i+1,4} = imgIDs{i}{3};
                        tempT{i+1,5} = imgIDs{i}{4};
                        tempT{i+1,6} = imgIDs{i}{2};
%                         tempT{i+1,7} = [imgIDs{i}{1} '_' imgIDs{i}{3} '_' imgIDs{i}{4}]; % use to link the DIC to the movies
                        tempT{i+1,7} = [imgIDs{i}{1} '_' imgIDs{i}{3}]; % use to link the DIC to the movies
                    end
                    % get the imaging period (T) and frequency (Fs) from the file
                    imgInfo = imfinfo(fullfile(imgFiles(imgFltr(i)).folder, imgFiles(imgFltr(i)).name));
                    switch app.options.Microscope
                        case 'Nikon A1'
                            T = imgInfo(1).ImageDescription;
                            T = regexp(T, '=', 'split');
                            T = sscanf(T{6}, '%f');
                            if isinf(T)
                                T = 1/(app.options.Frequency);
                            end
                        case 'Others'
                            if app.options.Frequency > 0
                                T = 1/(app.options.Frequency);
                            else
                                answerInfo = inputdlg({'Imaging frequecy (Hz)'}, 'Import data');
                                T = 1/str2double(answerInfo{1});
                            end
                    end
                    tempT{i+1,8} = [imgInfo(1).Width, imgInfo(1).Height, length(imgInfo), 1/T];
                    if strcmp(imgInfo(1).ByteOrder, 'big-endian')
                        tempT{i+1,9} = {[imgInfo.StripOffsets]+1};
                    else
                        tempT{i+1,9} = {[imgInfo.StripOffsets]};
                    end
                end
                delete(hWait);
                app.imgT = cell2table(tempT(2:end,:), 'VariableNames', tempT(1,:));
                % Select and show the first DIC image
                app.curDIC = app.dicT.CellID{1};
                detrendData(app);
                updateDIC(app, true)
                ListImagesChange(app)
                app.ToggleViewStack.Enable = 'on';
            else
                togglePointer(app)
                return
            end
            togglePointer(app)
        end
        
        function FileMenuLoadSelected(app, event)
            [fileName, filePath] = uigetfile(app.options.LastPath, 'Select Network File');
            if filePath ~= 0
                % Save the path to the settings
                app.options.LastPath = filePath;
                togglePointer(app)
                networkFiles = load(fullfile(filePath, fileName));
%                 if isfield(networkFiles, 'imgStore')
%                     if ~isfile(networkFiles.imgStore.Files(1))
%                         % Files are re-located
%                         newPath = uigetdir(app.options.LastPath, 'Relocate files');
%                         networkFiles.imgStore = imageDatastore(newPath, 'FileExtensions', {'.tif'});
%                         for ii = 1:size(networkFiles.dicT,1)
%                             [~,imgName,imgExt] = fileparts(networkFiles.dicT.Filename{ii});
%                             %networkFiles.dicT.Filename(ii) = fullfile(networkFiles.imgStore.Folders, [imgName, imgExt]);
%                             networkFiles.dicT.Filename(ii) = networkFiles.imgStore.Files(contains(networkFiles.imgStore.Files, imgName));
%                         end
%                         for ii = 1:size(networkFiles.imgT,1)
%                             [~,imgName,imgExt] = fileparts(networkFiles.imgT.Filename{ii});
%                             %networkFiles.imgT.Filename(ii) = fullfile(networkFiles.imgStore.Folders, [imgName, imgExt]);
%                             networkFiles.imgT.Filename(ii) = networkFiles.imgStore.Files(contains(networkFiles.imgStore.Files, imgName));
%                         end
%                         %dicFltr = contains(networkFiles.imgStore.Files, app.options.StillName);
%                         %networkFiles.dicT.Filename = networkFiles.imgStore.Files(dicFltr);
%                         %networkFiles.imgT.Filename = networkFiles.imgStore.Files(~dicFltr);
%                     end
%                     if isempty(app.imgDatastore)
%                         app.imgDatastore = networkFiles.imgStore;
%                     else
%                         oldImgPath = unique(cellfun(@fileparts, app.imgDatastore.Files, 'UniformOutput', false));
%                         allImgPaths = [oldImgPath; fileparts(networkFiles.imgStore.Files{1})];
%                         app.imgDatastore = imageDatastore(allImgPaths);
%                     end
%                 end
                if isfield(networkFiles, 'dicT')
                    if isempty(app.dicT)
                        app.dicT = networkFiles.dicT;
                    else
                        app.dicT = [app.dicT; networkFiles.dicT];
                    end
                    app.DicMenu.Enable = true;
                    app.bcDIC.Visible = true;
                end
                if isfield(networkFiles, 'imgT')
                    if isempty(app.imgT)
                        app.imgT = networkFiles.imgT;
                    else
                        if numel(app.imgT.Properties.VariableNames) ~= numel(networkFiles.imgT.Properties.VariableNames)
                            % We need to add columns
                            missVars = setdiff(app.imgT.Properties.VariableNames, networkFiles.imgT.Properties.VariableNames);
                            for missVar = missVars
                                networkFiles.imgT = addvars(networkFiles.imgT, nan(size(networkFiles.imgT,1), 1), 'NewVariableNames',missVar);
                            end
                        end
                        app.imgT = [app.imgT; networkFiles.imgT];
                    end
                    app.AnalysisMenu.Enable = true;
                    app.ToggleViewStack.Enable = 'on';
                end
                if isfield(networkFiles, 'fullT')
                    app.fullT = networkFiles.fullT;
                end
                if isfield(networkFiles, 'dicT') && isfield(networkFiles, 'imgT')
                    app.RadioSingleTrace.Enable = 'on';
                    app.RadioAllMean.Enable = 'on';
                    app.PlotMenu.Enable = 'on';
                    app.CellNumberText.Enable = 'on';
                    app.FileMenuExport.Enable = 'on';
                    app.curDIC = app.dicT.CellID{1};
                    if ~strcmp('DetrendData',app.imgT.Properties.VariableNames)
                        detrendData(app)
                    end
                    updateDIC(app, true)
                else
                    populateImageList(app)
                end
            else
                return
            end
            togglePointer(app)
        end
        
        function FileMenuSaveSelected(app, event)
            % First save the settings
            saveSettings(app)
            % Then save the data
            oldDir = cd(app.options.LastPath);
            [fileName, filePath] = uiputfile('*.mat', 'Save network data');
            savePath = fullfile(filePath, fileName);
            dicT = app.dicT;
            imgT = app.imgT;
%             imgStore = app.imgDatastore;
            fullT = app.fullT;
            save(savePath, 'dicT', 'imgT', 'fullT');
            cd(oldDir)
        end
        
        function FileMenuExportSelected(app, event)
            % Ask what needs to be saved
            whatExport = questdlg('What would you like to export?', 'Export as csv', 'Analysis', 'Traces', 'Both', 'Analysis');
            imgT = app.imgT;
            imgT.ImgProperties = imgT.ImgProperties(:,4);
            imgT.Properties.VariableNames{8} = 'Fs';
            imgT = imgT(:, [2:8 17 21 31:53]);
            switch whatExport
                case 'Analysis'
                    [fileName, filePath] = uiputfile('*.csv', 'Export network data');
                    writetable(imgT, fullfile(filePath, fileName));
                case 'Traces'
                    [fileName, filePath] = uiputfile('*.xlsx', 'Export network traces');
                    hWait = waitbar(0, 'Exporting data');
                    nSheet = size(app.imgT, 1);
                    allRaw = app.imgT.RawIntensity;
                    allFF0 = app.imgT.FF0Intensity;
                    allCellID = app.imgT.CellID;
                    %timeEst = 1;
                    for s = 1:nSheet
                        tic;
                        waitbar(s/nSheet, hWait, sprintf('Exporting data (~%.2f s)', timeEst*(nSheet-s)));
                        writematrix(allRaw{s}, fullfile(filePath, sprintf('Raw_%s', fileName)), 'Sheet', allCellID{s});
                        writematrix(allFF0{s}, fullfile(filePath, sprintf('FF0_%s', fileName)), 'Sheet', allCellID{s});
                        timeEst = toc;
                    end
                    close(hWait)
                otherwise
                    [fileName, filePath] = uiputfile('*.csv', 'Export network data');
                    writetable(imgT, fullfile(filePath, fileName));
                    fileName = regexprep(fileName, 'csv', 'xlsx');
                    hWait = waitbar(0, 'Exporting data');
                    nSheet = size(app.imgT, 1);
                    allRaw = app.imgT.RawIntensity;
                    allFF0 = app.imgT.FF0Intensity;
                    allCellID = app.imgT.CellID;
                    %timeEst = 1;
                    for s = 1:nSheet
                        tic;
                        waitbar(s/nSheet, hWait, sprintf('Exporting data (~%.2f s)', timeEst*(nSheet-s)));
                        writematrix(allRaw{s}, fullfile(filePath, sprintf('Raw_%s', fileName)), 'Sheet', allCellID{s});
                        writematrix(allFF0{s}, fullfile(filePath, sprintf('FF0_%s', fileName)), 'Sheet', allCellID{s});
                        timeEst = toc;
                    end
                    close(hWait)
            end
            
            
        end
        
        function DicMenuSelectSelected(app, event)
            % get a list of DIC images in the datastore
            dicFiles = app.dicT.CellID;
            [dicIndex, isOKed] = listdlg('PromptString', 'Please select DIC image',...
                'ListString', dicFiles, 'InitialValue', find(contains(dicFiles, app.curDIC)),...
                'SelectionMode', 'single');
            if isOKed
                dicName = dicFiles{dicIndex};
                app.curDIC = dicName;
            else
                return
            end
            bTable = false;
            if ~isempty(app.dicT) && sum(contains(app.dicT.CellID, dicName))
                bTable = true;
            end
            togglePointer(app)
            updateDIC(app, bTable)
            ListImagesChange(app)
            togglePointer(app)
        end
        
        function DicMenuPlaceRoiSelected(app, event)
            % Allow the user to place ROI on the DIC image
            axes(app.AxesDIC);
            [c, r, ~]= impixel;
            newRoi = [r, c];
            newRoi(newRoi <= app.options.RoiSize) = app.options.RoiSize + 1;
            newRoi(newRoi >= app.hDIC.XData(2) - app.options.RoiSize) = app.hDIC.XData(2) - app.options.RoiSize - 1;
            if ~isempty(app.dicT{contains(app.dicT.CellID, app.curDIC), 'RoiSet'})
                newRoi = [cell2mat(app.dicT{contains(app.dicT.CellID, app.curDIC), 'RoiSet'}); newRoi];
            end
            app.dicT(contains(app.dicT.CellID, app.curDIC), 'RoiSet') = {newRoi};
            % perform a basic analysis (raw intensity, FF0 data)
            updateDIC(app, true)
            getIntensityvalues(app)
            updatePlot(app)
            app.AnalysisMenu.Enable = true;
            app.RadioSingleTrace.Enable = 'on';
            app.RadioAllMean.Enable = 'on';
        end
        
        function DicMenuRemoveRoiSelected(app, event)
            % Allow the user to remove ROI on the DIC image
            axes(app.AxesDIC);
            [c, r, ~]= impixel;
            currRois = app.cellsMask;
            keepRoi = true(size(currRois,1),1);
            for cD = 1:numel(c)
                for roi = 1:size(currRois,1)
                    if inpolygon(r(cD),c(cD),app.cellsMask{roi,1}, app.cellsMask{roi,2})
                        keepRoi(roi) = false;
                    end
                end
            end
            currRois = cell2mat(app.dicT{contains(app.dicT.CellID, app.curDIC), 'RoiSet'});
            app.dicT{contains(app.dicT.CellID, app.curDIC), 'RoiSet'} = {currRois(keepRoi,:)};
            app.imgT{contains(app.imgT.CellID, app.curStak), 'RawIntensity'} = {app.imgT.RawIntensity{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
            app.imgT{contains(app.imgT.CellID, app.curStak), 'FF0Intensity'} = {app.imgT.FF0Intensity{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
            app.imgT{contains(app.imgT.CellID, app.curStak), 'DetrendData'} = {app.imgT.DetrendData{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
            if ~isempty(app.imgT{contains(app.imgT.CellID, app.curStak), 'SpikeLocations'})
                app.imgT{contains(app.imgT.CellID, app.curStak), 'SpikeLocations'} = {app.imgT.SpikeLocations{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
                app.imgT{contains(app.imgT.CellID, app.curStak), 'SpikeIntensities'} = {app.imgT.SpikeIntensities{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
                app.imgT{contains(app.imgT.CellID, app.curStak), 'SpikeWidths'} = {app.imgT.SpikeWidths{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
                app.imgT{contains(app.imgT.CellID, app.curStak), 'SpikeRaster'} = {app.imgT.SpikeRaster{contains(app.imgT.CellID, app.curStak)}(keepRoi,:)};
            end
            updateDIC(app, true)
            %getIntensityvalues(app)
            updatePlot(app)
        end
        
        function AnalysisMenuDetectSelected(app, event)
            answer = questdlg('Detect spike in which set?', 'Spike detection', 'All FOVs', 'Current set', 'Current FOV', 'All FOVs');
            switch answer
                case 'All FOVs'
                    imgFltr = find(~cellfun(@isempty, app.imgT{:,'RawIntensity'}));
                case 'Current set'
                    dicID = app.dicT.ExperimentID{contains(app.dicT.CellID, app.curDIC)};
                    imgFltr = find(strcmp(app.imgT.ExperimentID, dicID));
                case 'Current FOV'
                    imgFltr = find(strcmp(app.imgT.CellID, app.curStak));
            end
            nImg = numel(imgFltr);
            hWait = waitbar(0, 'Detecting spike in data');
            for i = 1:nImg
                waitbar(i/nImg, hWait, sprintf('Detecting spike in data %0.2f%%', i/nImg*100));
                detectSpike(app, imgFltr(i))
            end
            delete(hWait)
            app.TogglePeakRemove.Enable = 'on';
            app.TogglePeakAdd.Enable = 'on';
            app.RadioSingleTrace.Enable = 'on';
            app.RadioAllMean.Enable = 'on';
            app.CellNumberText.Enable = 'on';
        end
        
        function AnalysisMenuQuantifySelected(app, event)
            imgFltr = find(~cellfun(@isempty, app.imgT{:,'SpikeLocations'}));
            nImg = numel(imgFltr);
            hWait = waitbar(0, 'Collecting spike data');
            timeEst = 10 * nImg;
            for i = 1:nImg
                if timeEst > 60
                    waitbar(i/nImg, hWait, sprintf('Quantifying %s (~ %0.2f min)', regexprep(app.imgT.CellID{i}, '_', ' '), timeEst/60));
                else    
                    waitbar(i/nImg, hWait, sprintf('Quantifying %s (~ %0.2f s)', regexprep(app.imgT.CellID{i}, '_', ' '), timeEst));
                end
                tic
                quantifySpikes(app, imgFltr(i))
                timeEst = mean([timeEst, round(toc * (nImg - i))]);
            end
            delete(hWait)
            % Add the median per spike, consider that not all cells are quantified yet
            for isi = 1:size(app.imgT,1)
                if ~isempty(app.imgT.SpikeIntensities{isi})
                    % Suprathreshold
                    app.imgT.MeanInt(isi) = nanmean(cellfun(@(x,y) nanmean(x(~y)), app.imgT.SpikeIntensities{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanFWHM(isi) = nanmean(cellfun(@(x,y) nanmean(x(~y)), app.imgT.SpikeWidths{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanRiseTime(isi) = nanmean(cellfun(@(x,y) nanmean(x(~y)), app.imgT.CellRiseTime{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanDecayTau(isi) = nanmean(cellfun(@(x,y) nanmean(x(~y)), app.imgT.CellDecayTau{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianInt(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(~y)), app.imgT.SpikeIntensities{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianFWHM(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(~y)), app.imgT.SpikeWidths{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianRiseTime(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(~y)), app.imgT.CellRiseTime{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianDecayTau(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(~y)), app.imgT.CellDecayTau{isi}, app.imgT.SubthresholdFltr{isi}));
                    % Subthreshold
                    app.imgT.MeanSubInt(isi) = nanmean(cellfun(@(x,y) nanmean(x(y)), app.imgT.SpikeIntensities{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanSubFWHM(isi) = nanmean(cellfun(@(x,y) nanmean(x(y)), app.imgT.SpikeWidths{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanSubRiseTime(isi) = nanmean(cellfun(@(x,y) nanmean(x(y)), app.imgT.CellRiseTime{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MeanSubDecayTau(isi) = nanmean(cellfun(@(x,y) nanmean(x(y)), app.imgT.CellDecayTau{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianSubInt(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(y)), app.imgT.SpikeIntensities{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianSubFWHM(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(y)), app.imgT.SpikeWidths{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianSubRiseTime(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(y)), app.imgT.CellRiseTime{isi}, app.imgT.SubthresholdFltr{isi}));
                    app.imgT.MedianSubDecayTau(isi) = nanmedian(cellfun(@(x,y) nanmedian(x(y)), app.imgT.CellDecayTau{isi}, app.imgT.SubthresholdFltr{isi}));
                    % Interspike interval parameters
                    app.imgT.MeanISI(isi) = nanmean(cellfun(@nanmean, app.imgT.InterSpikeInterval{isi}));
                    app.imgT.MedianISI(isi) = nanmedian(cellfun(@nanmedian, app.imgT.InterSpikeInterval{isi}));
                end
            end
            app.imgT.ISICoV = cellfun(@nanmean, app.imgT.InterSpikeIntervalCoV);
            % For the Gaussian fit
            app.imgT.MeanFreq = cellfun(@nanmean, app.imgT.CellFrequency);
            app.imgT.MedianFreq = cellfun(@nanmedian, app.imgT.CellFrequency);
            app.imgT.nCell = cellfun(@(x) size(x,1), app.imgT.RawIntensity);
            app.imgT.nCell = round(app.imgT.nCell .* app.imgT.Participation / 100);
            app.imgT.PeakParticipation = cellfun(@nanmean, app.imgT.NetworkPeaks) ./ app.imgT.nCell * 100;
            app.PlotMenu.Enable = 'on';
            app.FileMenuExport.Enable = 'on';
        end
        
        function AnalysisMenuReCalculateSelected(app, event)
            doneFltr = find(~cellfun(@isempty, app.dicT{:,'RoiSet'}));
            startDIC = find(strcmp(app.dicT.CellID, app.curDIC));
            for i = 1:numel(doneFltr)
                app.curDIC = app.dicT{doneFltr(i), 'CellID'};
                getCellMask(app, doneFltr(i))
                getIntensityvalues(app)
                updatePlot(app)
            end
            app.curDIC = app.dicT{startDIC, 'CellID'};
            getCellMask(app, startDIC)
        end
        
        function PlotRasterSelected(app, event)
            bTrace = false;
            imgID = contains(app.imgT.CellID, app.curStak);
            Fs = app.imgT.ImgProperties(imgID,4);
            netRaster = app.imgT.NetworkRaster{imgID};
            nFrames = size(netRaster,2);
            time = (0:nFrames-1)/Fs;
            switch event.Source.Text
                case 'Event duration'
                    tempRaster = app.imgT.SpikeRaster{imgID};
                case 'Spike Locations'
                    tempLocs = app.imgT.SpikeLocations{imgID};
                    % Create the raster plot
                    nCell = size(tempLocs, 1);
                    tempRaster = nan(nCell, nFrames);
                    for c = 1:nCell
                        sStart = ceil(tempLocs{c} * Fs);
                        if ~isempty(sStart)
                            for s = 1:length(sStart)
                                tempRaster(c,sStart(s)-2:sStart(s)+2) = c;
                            end
                        end
                    end
                case 'Traces'
                    bTrace = true;
                    tempRaster = app.imgT.DetrendData{imgID};
                    cellSpace = max(tempRaster,[],'all') / 3;
                    cellNum = (1:size(tempRaster,1)) * cellSpace;
                    tempRaster = (tempRaster + repmat(cellNum',1,size(tempRaster,2)))';
            end
            % Actual plot
            figure()
            rAx = axes();
            hold on
            if bTrace
                plot(rAx, time, tempRaster, 'k')
                yMin = round(min(tempRaster,[],'all'), 2, 'significant');
                yMax = round(max(tempRaster,[],'all'), 2, 'significant');
                rAx.YLim = [yMin, yMax];
                rAx.YTick = linspace(yMin+cellSpace, yMax-cellSpace, size(tempRaster,2));
                rAx.YTickLabel = 1:size(tempRaster,2);
            else
                area(rAx, time, netRaster, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', .3)
                plot(rAx, time, tempRaster, 'k', 'LineWidth', 2)
                plot(rAx, [time(1) time(end)], ones(2,1) * sum(any(~isnan(tempRaster), 2)) * 0.8, 'r')
                rAx.YLim = [0, size(tempRaster,1)];
                rAx.YTick = 1:size(tempRaster,1);
            end
            rAx.Box = 'off';
            rAx.TickDir = 'out';
            rAx.YLabel.String = 'Cell number';
            rAx.XLabel.String = 'Time (s)';
            rAx.Title.String = regexprep(app.curStak, '_', ' ');
        end
        
        function OptionDebugSelected(app, event)
            disp(['You are now in debug mode. To exit the debug use "dbquit".',...
                ' Use "dbcont" to continue with the changes made'])
            keyboard
        end
        
        function ListImagesChange(app, event)
            warning('off', 'all')
            if numel(app.ListImages.String) < app.ListImages.Value
                app.ListImages.Value = 1;
            end
            app.curStak = app.ListImages.String{app.ListImages.Value};
            togglePointer(app)
            imgID = contains(app.imgT.Filename, app.curStak);
            tempData = app.imgT.DetrendData{imgID};
            yMin = min(tempData, [], 'all');
            yMax = max(tempData, [], 'all');
            app.YMinMax = [yMin, yMax];
            updatePlot(app)
            if app.ToggleViewStack.Value == 1
                createStackMovie(app);
            else
                togglePointer(app)
                imgFile = app.imgT.Filename{contains(app.imgT.Filename, app.curStak)};
                [~, fileName] = fileparts(imgFile);
                if isfile(imgFile)
                    imgTif = imread(imgFile);
                    fileFltr = find(contains(app.imgT.CellID, fileName));
                    if ~contains(fileName, app.options.Reference) && app.options.Registration
                        if any(strcmp('RegEst', app.imgT.Properties.VariableNames))
                            if isempty(cell2mat(app.imgT{fileFltr, 'RegEst'}))
                                tyrFile = app.imgT{find(contains(app.imgT{1:fileFltr,'CellID'}, 'Tyrode'),1,'last'),'Filename'};
                                imgTyr = imread(cell2mat(tyrFile));
                                imgOut = imref2d(size(imgTyr));
                                tformEstimate = imregcorr(imgTif, imgTyr, 'translation');
                                app.imgT{fileFltr,'RegEst'} = {tformEstimate.T};
                            else
                                tformEstimate = affine2d(cell2mat(app.imgT{fileFltr,'RegEst'}));
                                imgOut = imref2d(size(imgTif));
                            end
                        end
                        imgTif = imwarp(imgTif, tformEstimate, 'OutputView', imgOut);
                    end
                    app.minMax.Stack = stretchlim(imgTif,[0 .999]);
                    app.hStack = imshow(imadjust(imgTif, app.minMax.Stack), 'Parent', app.AxesStack);
                    app.movieData = [];
                    app.curSlice = [];
                    app.AxesStack.Title.String = regexprep(app.curStak, '_', ' ');
                else
                    app.AxesStack.Title.String = sprintf('File %s not found', regexprep(app.curStak, '_', ' '));
                    app.AxesStack.Visible = 'on';
                end
                togglePointer(app)
            end
            warning('on', 'all')
            % Add the mask to the movie
            if ~isempty(app.cellsMask)
                copyobj(app.patchMask, app.AxesStack)
            end
            togglePointer(app)
        end
        
        function ButtonPlotChange(app, event)
            switch event.NewValue.String
                case 'Single traces'
                    app.CellNumberMinor.Enable = 'on';
                    app.CellNumberText.Enable = 'on';
                    app.CellNumberPlus.Enable = 'on';
                    app.ToggleFixAxis.Enable = 'on';
                    if ~isempty(app.imgT.SpikeLocations{contains(app.imgT.CellID, app.curStak)})
                        app.TogglePeakAdd.Enable = 'on';
                        app.TogglePeakRemove.Enable = 'on';
                    end
                case 'All and mean'
                    app.CellNumberMinor.Enable = 'off';
                    app.CellNumberText.Enable = 'off';
                    app.CellNumberPlus.Enable = 'off';
                    app.TogglePeakAdd.Enable = 'off';
                    app.TogglePeakRemove.Enable = 'off';
                    app.ToggleFixAxis.Enable = 'off';
            end
            togglePointer(app)
            updatePlot(app)
            togglePointer(app)
        end
        
        function CellNumberEntered(app, event)
            oldCell = app.curCell(1);
            app.curCell = str2double(app.CellNumberText.String);
            app.curCell(2) = oldCell;
            togglePointer(app)
            updatePlot(app)
            togglePointer(app)
        end
        
        function CellNumberPlusPress(app, event)
            oldCell = app.curCell(1);
            nCells = numel(app.patchMask);
            if oldCell+1 > nCells
                app.CellNumberText.String = 1;
            else
                app.CellNumberText.String = oldCell+1;
            end
            app.curCell = str2double(app.CellNumberText.String);
            app.curCell(2) = oldCell;
            togglePointer(app)
            updatePlot(app)
            togglePointer(app)
        end
        
        function CellNumberMinorPress(app, event)
            oldCell = app.curCell(1);
            nCells = numel(app.patchMask);
            if oldCell-1 < 1
                app.CellNumberText.String = nCells;
            else
                app.CellNumberText.String = oldCell-1;
            end
            app.curCell = str2double(app.CellNumberText.String);
            app.curCell(2) = oldCell;
            togglePointer(app)
            updatePlot(app)
            togglePointer(app)
        end
        
        function SliderImageStackMoved(app, event)
            if app.SliderImageStack.Visible
                if nargin == 2 && isprop(event, 'VerticalScrollCount')
                    sliceToShow = round(app.SliderImageStack.Value + event.VerticalScrollCount);
                    if sliceToShow < 1
                        return
                    end
                else
                    sliceToShow = round(app.SliderImageStack.Value);
                end
                app.SliderImageStack.Value = sliceToShow;
                Fs = app.imgT.ImgProperties(contains(app.imgT.CellID, app.curStak),4);
                app.curSlice.CData = app.movieData(:,:,sliceToShow);
                app.curTime.XData = ones(2,1)*sliceToShow/Fs - 1/Fs;
            end
        end
        
        function ToggleViewStackPressed(app, event)
            if app.ToggleViewStack.Value == 1
                if isempty(app.movieData) || ~matches(app.AxesStack.Title.String, regexprep(app.curStak, '_', ' '))
                    if ~isfile(app.imgT.Filename(1))
                        % Files are re-located
                        newPath = uigetdir(app.options.LastPath, 'Relocate files');
                        dicFiles = fullfile(newPath, {app.dicT.CellID{:}})';
                        app.dicT.Filename = cellfun(@(x) sprintf('%s.tif', x), dicFiles, 'UniformOutput', false);
                        imgFiles = fullfile(newPath, {app.imgT.CellID{:}})';
                        app.imgT.Filename = cellfun(@(x) sprintf('%s.tif', x), imgFiles, 'UniformOutput', false);
                    end
                    togglePointer(app)
                    app.ToggleViewStack.BackgroundColor = 'white';
                    createStackMovie(app)
                    togglePointer(app)
                end
            else
                app.ToggleViewStack.BackgroundColor = [.94 .94 .94];
                app.SliderImageStack.Visible = 'off';
                app.movieData = [];
                app.curSlice = [];
                app.curTime.Visible = 'off';
                app.curTime = [];
            end
        end
        
        function TogglePeakRemovePress(app, event)
            % First get the filters for the data
            imgFltr = find(contains(app.imgT.CellID, app.curStak));
            cellFltr = app.curCell(1);
            Fs = app.imgT.ImgProperties(imgFltr,4);
            spikeLoc = app.imgT{imgFltr, 'SpikeLocations'}{:};
            spikeInt = app.imgT{imgFltr, 'SpikeIntensities'}{:};
            spikeWid = app.imgT{imgFltr, 'SpikeWidths'}{:};
            spikeData = spikeLoc{cellFltr,1};
            keepEvent = true(size(spikeData));
            while true
                axes(app.AxesPlot);
                [pointX, ~] = ginput(1);
                sel = get(app.Figure, 'SelectionType');
                if strcmpi(sel, 'alt')
                    break
                end
                %pointX = round(round(pointX) * Fs);
                delPeak = find(spikeData > pointX-30/Fs & spikeData < pointX+30/Fs);
                if numel(delPeak) > 1
                    delPeak = find(spikeData > pointX-10/Fs & spikeData < pointX+10/Fs);
                    if numel(delPeak) > 1
                        % Narrow down to 2 frames
                        delPeak = find(spikeData > pointX-2/Fs & spikeData < pointX+2/Fs);
                        if numel(delPeak) > 1
                            % Get the first one in the list
                            delPeak = delPeak(1);
                        end
                    end
                end
                if numel(delPeak) > 1
                    hWarn = warndlg('Multiple peaks detected at this location', 'Delete peak');
                    waitfor(hWarn);
                elseif numel(delPeak) == 0
                    hWarn = warndlg('No peaks detected at this location', 'Delete peak');
                    waitfor(hWarn);
                else
                    if keepEvent(delPeak) == 0
                        keepEvent(delPeak) = 1;
                        currDel = findobj(app.AxesPlot, 'LineWidth', 1.5, 'XData', spikeLoc{cellFltr, 1}(delPeak));
                        currDel.Visible = 'off';
                    else
                        keepEvent(delPeak) = 0;
                        currDel = findobj(app.AxesPlot, 'LineWidth', 1.5, 'XData', spikeLoc{cellFltr, 1}(delPeak));
                        if ~isempty(currDel)
                            currDel.Visible = 'on';
                        end
                        plot(app.AxesPlot, spikeLoc{cellFltr, 1}(delPeak), spikeInt{cellFltr, 1}(delPeak), 'xr', 'LineWidth', 1.5)
                    end
                end
            end
            % Adjust the data and show the new plot
            app.imgT{imgFltr, 'SpikeLocations'}{1}{cellFltr} = spikeLoc{cellFltr, 1}(keepEvent);
            app.imgT{imgFltr, 'SpikeIntensities'}{1}{cellFltr} = spikeInt{cellFltr, 1}(keepEvent);
            app.imgT{imgFltr, 'SpikeWidths'}{1}{cellFltr} = spikeWid{cellFltr, 1}(keepEvent);
            updatePlot(app);
        end
        
        function TogglePeakAddPress(app, event)
            % First get the filters for the data
            imgFltr = find(contains(app.imgT.CellID, app.curStak));
            cellFltr = app.curCell(1);
            Fs = app.imgT.ImgProperties(imgFltr,4);
            spikeLoc = app.imgT{imgFltr, 'SpikeLocations'}{:};
            spikeInt = app.imgT{imgFltr, 'SpikeIntensities'}{:};
            spikeWid = app.imgT{imgFltr, 'SpikeWidths'}{:};
            spikeData(1,:) = spikeLoc{cellFltr,1};
            spikeData(2,:) = spikeInt{cellFltr,1};
            spikeData(3,:) = spikeWid{cellFltr,1};
            tempSpikeData = nan(3,1);
            % Get the trace and the detection paramethers
            trace = cell2mat(app.imgT{imgFltr, 'DetrendData'});
            trace = trace(cellFltr,:);
            smoothTrace = wdenoise(trace, 'DenoisingMethod', 'BlockJS');
            halfDuration = round((app.options.PeakMinDuration + app.options.PeakMaxDuration) / 2);
            spikeLeng = halfDuration / Fs;
            switch app.options.DetectTrace
                case 'Raw'
                    detectTrace = trace;
                case 'Gradient'
                    detectTrace = gradient(smoothTrace);
                case 'Smooth'
                    detectTrace = smoothTrace;
            end
            tS = 1;
            while true
                axes(app.AxesPlot);
                [pointX, ~] = ginput(1);
                sel = get(app.Figure, 'SelectionType');
                if strcmp(event.EventName, 'KeyPress')
                    currDel = findobj(app.AxesPlot, 'LineWidth', 1.5, 'XData', tempSpikeData(2,tS));
                    tempSpikeData(1,tS-1) = [];
                    tempSpikeData(2,tS-1) = [];
                    tempSpikeData(3,tS-1) = [];
                else
                    if strcmpi(sel, 'alt')
                        break
                    end
                    if isempty(pointX)
                        tempSpikeData = nan;
                        break
                    end
                    pointF = round(round(pointX) * Fs);
                    [tempInts, tempLocs, tempWidths] = findpeaks(detectTrace(pointF-halfDuration:pointF+halfDuration),...
                        'SortStr', 'descend');
                    tempSpikeData(1,tS) = (pointF-halfDuration+tempLocs(1)-2) / Fs;
                    % Check to see that the new spike is at least as far as the minimum distance
                    fltr1 = (tempSpikeData(1,tS)-app.options.PeakMinDistance) < spikeData(1,:);
                    fltr2 = spikeData(1, fltr1) <= tempSpikeData(1,tS);
                    if any(fltr2)
                        hErr = errordlg('Spike already selected!');
                        waitfor(hErr);
                        tempSpikeData(:,tS) = [];
                    else
                        if any(tempSpikeData(1,1:tS-1)==tempSpikeData(1,tS))
                            hErr = errordlg('Spike already selected!');
                            waitfor(hErr);
                            tempSpikeData(:,tS) = [];
                        else
                            tempSpikeData(2,tS) = tempInts(1);
                            tempSpikeData(3,tS) = tempWidths(1);
                            plot(app.AxesPlot, tempSpikeData(1,tS), tempSpikeData(2,tS), 'og', 'LineWidth', 1.5)
                            tS = tS + 1;
                        end
                    end
                end
            end
            % Add the data to the table and sort in order of event
            if ~isnan(tempSpikeData(1))
                spikeData = [spikeData, tempSpikeData];
                [sortedData, sortIdx] = sort(spikeData(1,:));
                sortedData(2,:) = spikeData(2,sortIdx);
                sortedData(3,:) = spikeData(3,sortIdx);
                app.imgT{imgFltr, 'SpikeLocations'}{1}{cellFltr} = sortedData(1,:);
                app.imgT{imgFltr, 'SpikeIntensities'}{1}{cellFltr} = sortedData(2,:);
                app.imgT{imgFltr, 'SpikeWidths'}{1}{cellFltr} = sortedData(3,:);
            end
            updatePlot(app);
        end
        
        function ToggleFixAxisPress(app, event)
            if isempty(app.YMinMax)
                imgID = contains(app.imgT.Filename, app.curStak);
                tempData = app.imgT.DetrendData{imgID};
                yMin = min(tempData, [], 'all');
                yMax = max(tempData, [], 'all');
                app.YMinMax = [yMin, yMax];
            end
            if app.ToggleFixAxis.Value
                app.AxesPlot.YLim = app.YMinMax;
                app.ToggleFixAxis.BackgroundColor = [0 .47 .84];
            else
                app.AxesPlot.YLimMode = 'auto';
                app.ToggleFixAxis.BackgroundColor = [.94 .94 .94];
            end
                
        end
        
        function keyPressed(app, event)
            switch event.Key
                case "a" % Add new peaks
                    events.EventName='Action';
                    TogglePeakAddPress(app, events)
                case "s" % Select new DIC
                    DicMenuSelectSelected(app, event)
                case "d" % delete all events
                    bDel = questdlg('Delete all the event?', 'Delete events');
                    if strcmp(bDel, 'Yes')
                        imgFltr = find(contains(app.imgT.CellID, app.curStak));
                        cellFltr = app.curCell(1);
                        Fs = app.imgT.ImgProperties(imgFltr,4);
                        spikeLoc = app.imgT{imgFltr, 'SpikeLocations'}{:};
                        spikeInt = app.imgT{imgFltr, 'SpikeIntensities'}{:};
                        spikeWid = app.imgT{imgFltr, 'SpikeWidths'}{:};
                        spikeData = spikeLoc{cellFltr,1};
                        keepEvent = false(size(spikeData));
                        app.imgT{imgFltr, 'SpikeLocations'}{1}{cellFltr} = spikeLoc{cellFltr, 1}(keepEvent);
                        app.imgT{imgFltr, 'SpikeIntensities'}{1}{cellFltr} = spikeInt{cellFltr, 1}(keepEvent);
                        app.imgT{imgFltr, 'SpikeWidths'}{1}{cellFltr} = spikeWid{cellFltr, 1}(keepEvent);
                        updatePlot(app);
                    end 
                case "rightarrow" % move to next cell
                    if app.RadioSingleTrace.Value
                        CellNumberPlusPress(app, event)
                    end
                case "leftarrow" % move to previous cell
                    if app.RadioSingleTrace.Value
                        CellNumberMinorPress(app, event)
                    end
                case "uparrow" % move to previous FOV
                    if app.ListImages.Value > 1
                        app.ListImages.Value = app.ListImages.Value - 1;
                        ListImagesChange(app, event);
                    end
                case "downarrow" % move to next FOV
                    if app.ListImages.Value < numel(app.ListImages.String)
                        app.ListImages.Value = app.ListImages.Value + 1;
                        ListImagesChange(app, event);
                    end
                case "p" % place ROI
                    DicMenuPlaceRoiSelected(app, event)
                case "r" % Remove ROI
                    DicMenuRemoveRoiSelected(app, event)
            end
        end
        
    end
    
    % Create the UIFigure and components
    methods (Access = private)
        function createComponents(app)
            screenSize = get(0, 'ScreenSize');
            app.Figure = figure('Units', 'pixels', 'Visible', 'off',...
                'Position', [screenSize(3)/5 screenSize(4)/6 1250 850],...
                'Name', 'Network Activity', 'ToolBar', 'none', 'MenuBar', 'none',...
                'NumberTitle', 'off', 'WindowScrollWheelFcn', @(~,event)SliderImageStackMoved(app, event),...
                'KeyPressFcn', @(~,event)keyPressed(app, event));
            movegui(app.Figure, 'center');
            % Create the menu bar: file options
            app.FileMenu = uimenu(app.Figure, 'Text', 'File');
            app.FileMenuOpen = uimenu(app.FileMenu, 'Text', 'Load new data',...
                'MenuSelectedFcn', createCallbackFcn(app, @FileMenuOpenSelected, true));
            app.FileMenuLoad = uimenu(app.FileMenu, 'Text', 'Open saved data',...
                'MenuSelectedFcn', createCallbackFcn(app, @FileMenuLoadSelected, true));
            app.FileMenuSave = uimenu(app.FileMenu, 'Text', 'Save',...
                'MenuSelectedFcn', createCallbackFcn(app, @FileMenuSaveSelected, true));
            app.FileMenuLabelCondition = uimenu(app.FileMenu, 'Text', 'Label condition',...
                'MenuSelectedFcn', createCallbackFcn(app, @FileLabelConditionSelected, true));
            app.FileMenuExport = uimenu(app.FileMenu, 'Text', 'Export',...
                'MenuSelectedFcn', createCallbackFcn(app, @FileMenuExportSelected, true),...
                'Enable', 'off');
            % DIC menu
            app.DicMenu = uimenu(app.Figure, 'Text', 'DIC Menu', 'Enable', 'off');
            app.DicMenuPlaceRoi = uimenu(app.DicMenu, 'Text', 'Place ROIs',...
                'MenuSelectedFcn', createCallbackFcn(app, @DicMenuPlaceRoiSelected, true));
            app.DicMenuRemoveRoi = uimenu(app.DicMenu, 'Text', 'Remove ROIs',...
                'MenuSelectedFcn', createCallbackFcn(app, @DicMenuRemoveRoiSelected, true));
            app.DicMenuSelect = uimenu(app.DicMenu, 'Text', 'Select Image', 'Separator', 'on',...
                'MenuSelectedFcn', createCallbackFcn(app, @DicMenuSelectSelected, true));
            % Analysis menu
            app.AnalysisMenu = uimenu(app.Figure, 'Text', 'Analysis', 'Enable', 'off');
            app.AnalysisMenuDetect = uimenu(app.AnalysisMenu, 'Text', 'Detect Events',...
                'MenuSelectedFcn', createCallbackFcn(app, @AnalysisMenuDetectSelected, true));
            app.AnalysisMenuQuantify = uimenu(app.AnalysisMenu, 'Text', 'Quantify Events',...
                'MenuSelectedFcn', createCallbackFcn(app, @AnalysisMenuQuantifySelected, true));
            app.AnalysisMenuReCalculate = uimenu(app.AnalysisMenu, 'Text', 'Re-calculate traces',...
                'Separator', 'on', 'MenuSelectedFcn', createCallbackFcn(app, @AnalysisMenuReCalculateSelected, true));
            % Plot menu
            app.PlotMenu = uimenu(app.Figure, 'Text', 'Plot', 'Enable', 'off');
            app.PlotMenuRaster = uimenu(app.PlotMenu, 'Text', 'Raster Plot', 'Enable', 'on');
            app.PlotMenuRasterAll = uimenu(app.PlotMenuRaster, 'Text', 'Event duration',...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotRasterSelected, true), 'Enable', 'on');
            app.PlotMenuRasterSpike = uimenu(app.PlotMenuRaster, 'Text', 'Spike Locations',...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotRasterSelected, true), 'Enable', 'on');
            app.PlotMenuRasterTrace = uimenu(app.PlotMenuRaster, 'Text', 'Traces',...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotRasterSelected, true), 'Enable', 'on');
            app.PlotMenuOverview = uimenu(app.PlotMenu, 'Text', 'Plot overview',...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotOverviewSelected, true), 'Enable', 'on');
            app.PlotMenuTimeFrequency = uimenu(app.PlotMenu, 'Text', 'ISI over time',...
                'MenuSelectedFcn', createCallbackFcn(app, @PlotTimeFrequencySelected, true), 'Enable', 'on');
            % Option menu
            app.OptionMenu = uimenu(app.Figure, 'Text', 'Options');
            app.OptionMenuSettings = uimenu(app.OptionMenu, 'Text', 'Settings',...
                'MenuSelectedFcn', createCallbackFcn(app, @OptionMenuSelected, true), 'Enable', 'on');
            app.OptionMenuDebug= uimenu(app.OptionMenu, 'Text', 'Debug',...
                'MenuSelectedFcn', createCallbackFcn(app, @OptionDebugSelected, true), 'Enable', 'on');
            % Create the axes to store the images
            app.AxesDIC = axes('Units', 'pixels', 'Position', [70 410 400 400],...
                'Visible', 'off');
            app.AxesDIC.Title.String = 'DIC Image';
            app.bcDIC = uicontrol(app.Figure, 'Style', 'pushbutton', 'Units', 'pixels',...
                'Position', [70 810 80 20], 'String', 'Adjust B&C',...
                'Visible', 'off', 'Callback', createCallbackFcn(app, @bcDICpressed, true));
            app.AxesStack = axes('Units', 'pixels', 'Position', [550 410 400 400],...
                'Visible', 'off');
            app.manualReg = uicontrol(app.Figure, 'Style', 'pushbutton', 'Units', 'pixels',...
                'Position', [550 810 80 20], 'String', 'Registration',...
                'Visible', 'off', 'Callback', createCallbackFcn(app, @manualRegistrationPressed, true));
            app.AxesStack.Title.String = 'Ca^{2+} Movie';
            app.SliderImageStack = uicontrol(app.Figure, 'Style', 'slider',...
                'Units', 'pixels', 'Position', [550 385 400 20],...
                'Visible', 'off', 'Value', 1, 'SliderStep', [1 10],...
                'Callback', createCallbackFcn(app, @SliderImageStackMoved, true));
            addlistener(app.SliderImageStack,'Value','PreSet',@(~,~)SliderImageStackMoved(app));
            % Create the axes for the plots
            app.AxesPlot = axes('Units', 'pixels', 'Position', [70 55 890 300], 'Visible', 'on',...
                'ButtonDownFcn', createCallbackFcn(app, @PlotClicked, true));
            app.AxesPlot.Title.String = 'Ca^{2+} Traces';
            app.AxesPlot.YLabel.String = '\DeltaF/F_0';
            app.AxesPlot.XLabel.String = 'Time (s)';
            app.AxesPlot.TickDir = 'out';
            % Create the components to interact with the images
            app.ListImagesLabel = uicontrol(app.Figure, 'Style', 'text',...
                'String', 'Image Files', 'Units', 'pixels',...
                'Position', [1000 810 100 22], 'FontSize', 12);
            app.ListImages = uicontrol(app.Figure, 'Style', 'listbox',...
                'Units', 'pixels', 'Position', [1000 410 200 400],...
                'Callback', createCallbackFcn(app, @ListImagesChange));
            app.ToggleViewStack = uicontrol(app.Figure, 'Style', 'togglebutton',...
                'Units', 'pixels', 'Position', [1100 810 100 20],...
                'String', 'Show movie', 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @ToggleViewStackPressed, true));
            app.ButtonPlot = uibuttongroup(app.Figure, 'Units', 'pixels',...
                'Title', 'Plot Type', 'Position', [1000 280 170 80], 'FontSize', 12,...
                'SelectionChangedFcn', createCallbackFcn(app, @ButtonPlotChange, true));
            app.RadioAllMean = uicontrol(app.ButtonPlot, 'Style', 'radiobutton',...
                'String', 'All and mean', 'Units', 'pixels',...
                'Position', [10 30 110 30], 'FontSize', 12, 'Enable', 'off');
            app.RadioSingleTrace = uicontrol(app.ButtonPlot, 'Style', 'radiobutton',...
                'String', 'Single traces', 'Units', 'pixels',...
                'Position', [10 5 110 30], 'FontSize', 12, 'Enable', 'off');
            app.CellNumberLabel = uicontrol(app.Figure, 'Style', 'text',...
                'String', 'Cell number','Units', 'pixels',...
                'Position', [1000 235 110 30], 'FontSize', 12);
            app.CellNumberText = uicontrol(app.Figure, 'Style', 'edit',...
                'String', '1', 'Units', 'pixels', 'Position', [1050 215 50 30],...
                'FontSize', 12, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @CellNumberEntered, true));
            app.CellNumberPlus = uicontrol(app.Figure, 'Style', 'pushbutton',...
                'String', 'Next', 'Units', 'pixels', 'Position', [1100 215 50 30],...
                'FontSize', 12, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @CellNumberPlusPress, true));
            app.CellNumberMinor = uicontrol(app.Figure, 'Style', 'pushbutton',...
                'String', 'Prev', 'Units', 'pixels', 'Position', [1000 215 50 30],...
                'FontSize', 12, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @CellNumberMinorPress, true));
            app.TogglePeakAdd = uicontrol(app.Figure, 'Style', 'pushbutton',...
                'String', 'Add peaks', 'Units', 'pixels', 'Position', [1000 175 110 30],...
                'FontSize', 12, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @TogglePeakAddPress, true));
            app.TogglePeakRemove = uicontrol(app.Figure, 'Style', 'pushbutton',...
                'String', 'Remove peaks', 'Units', 'pixels', 'Position', [1000 135 110 30],...
                'FontSize', 12, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @TogglePeakRemovePress, true));
            app.ToggleFixAxis = uicontrol(app.Figure, 'Style', 'togglebutton',...
                'String', 'Fix Y Axis', 'Units', 'pixels', 'Position', [1000 95 110 30],...
                'FontSize', 13, 'Enable', 'off',...
                'Callback', createCallbackFcn(app, @ToggleFixAxisPress, true));
            % Set the figure to visible
            nUI = fieldnames(app);
            for ui = 1:numel(nUI)
                if ~contains(nUI{ui}, 'Menu')
                    app.(nUI{ui}).Units = 'normalized';
                end
            end
            app.Figure.Visible = 'on';
        end
        
        function optFigure = OptionMenuSelected(app, event)
            % Create a small GUI to select the options
            micOptions = {'Nikon A1', 'Others'};
            shapeOptions = {'Square', 'Circle'};
            heightOptions = {'MAD', 'Normalized MAD', 'Rolling StDev'};
            traceOptions = {'Raw', 'Gradient', 'Smooth'};
            detrendOptions = {'None', 'Moving median', 'Polynomial'};
            screenSize = get(0, 'ScreenSize');
            optFigure = figure('Units', 'pixels', 'Visible', 'on',...
                'Position', [screenSize(3)/2 screenSize(4)/2 500 500],...
                'Name', 'Network Activity Options', 'ToolBar', 'none', 'MenuBar', 'none',...
                'NumberTitle', 'off');
            pLoad = uipanel('Title', 'Load options', 'Position', [.01 .7 .4 .28], 'FontSize', 12);
                uicontrol('Parent', pLoad, 'Style', 'text', 'String', 'Still condition',...
                    'Position', [10 80 80 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hLoadStill = uicontrol('Parent', pLoad, 'Style', 'edit', 'String', app.options.StillName,...
                    'Position', [100 80 80 20], 'FontSize', 10);
                uicontrol('Parent', pLoad, 'Style', 'text', 'String', 'Microscope',...
                    'Position', [10 50 80 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hLoadMic = uicontrol('Parent', pLoad, 'Style', 'popupmenu', 'String', micOptions,...
                    'Value', find(contains(micOptions, app.options.Microscope)),...
                    'Position', [100 50 80 20], 'FontSize', 10);
                uicontrol('Parent', pLoad, 'Style', 'text', 'String', 'Frequency',...
                    'Position', [10 20 80 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hLoadFreq = uicontrol('Parent', pLoad, 'Style', 'edit', 'String', app.options.Frequency,...
                    'Position', [100 20 80 20], 'FontSize', 10);
            pRoi = uipanel('Title', 'ROIs options', 'Position', [.01 .4 .4 .28], 'FontSize', 12);
                uicontrol('Parent', pRoi, 'Style', 'text', 'String', 'ROI size (px)',...
                    'Position', [10 80 90 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hRoiSize = uicontrol('Parent', pRoi, 'Style', 'edit', 'String', app.options.RoiSize,...
                    'Position', [100 80 80 20], 'FontSize', 10);
                uicontrol('Parent', pRoi, 'Style', 'text', 'String', 'ROI shape',...
                    'Position', [10 50 80 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hRoiShape = uicontrol('Parent', pRoi, 'Style', 'popupmenu', 'String', shapeOptions,...
                    'Value', find(contains(shapeOptions, app.options.RoiShape)),...
                    'Position', [100 50 80 20], 'FontSize', 10);
            pReg = uipanel('Title', 'Registration options', 'Position', [.01 .1 .4 .28], 'FontSize', 12);
                hReg = uicontrol('Parent', pReg, 'Style', 'checkbox', 'String', 'Registration',...
                    'Value', app.options.Registration, 'Position', [5 80 100 20], 'FontSize', 10);
                uicontrol('Parent', pReg, 'Style', 'text', 'String', 'Reference condition',...
                    'Position', [5 50 150 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hRegRef = uicontrol('Parent', pReg, 'Style', 'edit', 'String', app.options.Reference,...
                    'Position', [5 30 150 20], 'FontSize', 10);
            pPeak = uipanel('Title', 'Detection options', 'Position', [.45 .38 .5 .6], 'FontSize', 12);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Height methods',...
                    'Position', [10 240 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hMinHeigth = uicontrol('Parent', pPeak, 'Style', 'popupmenu', 'String', heightOptions,...
                    'Value', find(contains(heightOptions, app.options.PeakMinHeight),1),...
                    'Position', [150 240 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Threshold',...
                    'Position', [10 210 100 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hSigma = uicontrol('Parent', pPeak, 'Style', 'edit', 'String', app.options.SigmaThr,...
                    'Position', [150 210 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Min Prominance',...
                    'Position', [10 180 100 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hMinProm = uicontrol('Parent', pPeak, 'Style', 'edit', 'String', app.options.PeakMinProminance,...
                    'Position', [150 180 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Min Distance (frames)',...
                    'Position', [10 150 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hMinDist = uicontrol('Parent', pPeak, 'Style', 'edit', 'String', app.options.PeakMinDistance,...
                    'Position', [150 150 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Min Duration (frames)',...
                    'Position', [10 120 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hMinDura = uicontrol('Parent', pPeak, 'Style', 'edit', 'String', app.options.PeakMinDuration,...
                    'Position', [150 120 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Max Duration (frames)',...
                    'Position', [10 90 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hMaxDura = uicontrol('Parent', pPeak, 'Style', 'edit', 'String', app.options.PeakMaxDuration,...
                    'Position', [150 90 80 20], 'FontSize', 10);
                uicontrol('Parent', pPeak, 'Style', 'text', 'String', 'Trace to use',...
                    'Position', [10 60 80 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hDetTrace = uicontrol('Parent', pPeak, 'Style', 'popupmenu', 'String', traceOptions,...
                    'Value', find(contains(traceOptions, app.options.DetectTrace),1),...
                    'Position', [150 60 80 20], 'FontSize', 10);
            pDetrend = uipanel('Title', 'Detrending', 'Position', [.45 .1 .5 .28], 'FontSize', 12);
                uicontrol('Parent', pDetrend, 'Style', 'text', 'String', 'Detrending method',...
                    'Position', [10 80 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hDetMeth = uicontrol('Parent', pDetrend, 'Style', 'popupmenu', 'String', detrendOptions,...
                    'Value', find(contains(detrendOptions, app.options.Detrending),1),...
                    'Position', [150 80 80 20], 'FontSize', 10);
                uicontrol('Parent', pDetrend, 'Style', 'text', 'String', 'Detrend window (fr)',...
                    'Position', [10 50 130 20], 'FontSize', 10, 'HorizontalAlignment', 'left');
                hDetWin = uicontrol('Parent', pDetrend, 'Style', 'edit', 'String', app.options.DetrendSize,...
                    'Position', [150 50 80 20], 'FontSize', 10);
            hSave = uicontrol(optFigure, 'Style', 'pushbutton',...
                'String', 'Save', 'Units', 'pixels', 'Position', [110 10 90 30],...
                'FontSize', 12, 'Enable', 'on',...
                'Callback', createCallbackFcn(app, @saveOptionPressed, true));
            hCancel = uicontrol(optFigure, 'Style', 'pushbutton',...
                'String', 'Cancel', 'Units', 'pixels', 'Position', [210 10 90 30],...
                'FontSize', 12, 'Enable', 'on',...
                'Callback', createCallbackFcn(app, @cancelOptionPressed, true));
            hReset = uicontrol(optFigure, 'Style', 'pushbutton',...
                'String', 'Reset', 'Units', 'pixels', 'Position', [310 10 90 30],...
                'FontSize', 12, 'Enable', 'on',...
                'Callback', createCallbackFcn(app, @resetOptionPressed, true));
            function cancelOptionPressed(app, event)
                delete(optFigure);
            end
            function saveOptionPressed(app, event)
                app.options.StillName = hLoadStill.String;
                app.options.Microscope = hLoadMic.String{hLoadMic.Value};
                app.options.Frequency = str2double(hLoadFreq.String);
                app.options.RoiSize = str2double(hRoiSize.String);
                app.options.RoiShape = hRoiShape.String{hRoiShape.Value};
                app.options.PeakMinHeight = hMinHeigth.String{hMinHeigth.Value};
                app.options.SigmaThr = str2double(hSigma.String);
                app.options.PeakMinProminance = str2double(hMinProm.String);
                app.options.PeakMinDistance = str2double(hMinDist.String);
                app.options.PeakMinDuration = str2double(hMinDura.String);
                app.options.PeakMaxDuration = str2double(hMaxDura.String);
                app.options.DetectTrace = hDetTrace.String{hDetTrace.Value};
                app.options.Registration = hReg.Value;
                app.options.Reference = hRegRef.String;
                app.options.Detrending = hDetMeth.String{hDetMeth.Value};
                app.options.DetrendSize = str2double(hDetWin.String);
                cancelOptionPressed(app, event);
                if ~isempty(app.imgT)
                    detrendData(app)
                    updatePlot(app);
                end
                if app.options.Registration
                    app.manualReg.Visible = 'on';
                else
                    app.manualReg.Visible = 'off';
                end
            end
            function resetOptionPressed(app, event)
                app.options.StillName = 'DIC.';
                app.options.Microscope = 'Nikon A1';
                app.options.Frequency = 0;
                app.options.RoiSize = 5;
                app.options.RoiShape = 'Square';
                app.options.PeakMinHeight = 'MAD';
                app.options.PeakMinProminance = 0.1;
                app.options.PeakMinDistance = 1;
                app.options.PeakMinDuration = 2;
                app.options.PeakMaxDuration = 5;
                app.options.DetectTrace = 'Smooth';
                app.options.Registration = false;
                app.options.Reference = 'Tyrode';
                app.options.Detrending = 'None';
                app.options.DetrendSize = 0;
                cancelOptionPressed(app, event);
                OptionMenuSelected(app, event);
            end
        end
        
        % Rename the conditions
        function FileLabelConditionSelected(app, ~)
            screenSize = get(0, 'ScreenSize');
            labelFigure = figure('Units', 'pixels', 'Visible', 'on',...
                'Position', [100 100 400 500],...
                'Name', 'Relable condition', 'ToolBar', 'none', 'MenuBar', 'none',...
                'NumberTitle', 'off');
            hList = uicontrol('Style', 'listbox', 'Position', [20 80 200 400],...
                'Min', 0, 'Max', 2, 'String', cellfun(@(x,y) strjoin({x,y}, ' - '), app.imgT.CellID, string(app.imgT.Condition), 'UniformOutput', false));
            hRadioButton = uibuttongroup('Title', 'Label to change',...
                'Units', 'pixels', 'Position', [230 409 160 80],...
                'FontSize', 11);
            hRadioCond = uicontrol(hRadioButton, 'Style', 'radiobutton',...
                'String', 'Condition', 'Position', [10 30 110 30]);
            hRadioRec = uicontrol(hRadioButton, 'Style', 'radiobutton',...
                'String', 'Recording', 'Position', [10 5 110 30]);
            hText = uicontrol('Style', 'edit', 'Position', [240 360 140 20],...
                'String', '(Condition)');
            hButton = uicontrol('Style', 'pushbutton', 'String', 'Confirm',...
                'Position', [240 335 140 20],...
                'Callback', createCallbackFcn(app, @labelConfirmPressed, true));
            function labelConfirmPressed(app, ~)
                % Get the index of the files to change
                nameIdx = hList.Value;
                newLabel = hText.String;
                switch hRadioButton.SelectedObject.String
                    case 'Condition'
                        app.imgT.Condition(nameIdx) = {newLabel};
                    case 'Recording'
                        app.imgT.RecID(nameIdx) = {newLabel};
                end
                hList.String = cellfun(@(x,y) strjoin({x,y}, ' - '), app.imgT.CellID, string(app.imgT.Condition), 'UniformOutput', false);
            end
            waitfor(labelFigure);
            % ask the user to order the conditions
            uniCond = string(unique(app.imgT.Condition));
            promptCond = 'Enter space-separated number for conditions: ';
            for c=1:numel(uniCond)
                promptCond = sprintf('%s, %s', promptCond, uniCond{c});
            end
            condOrd = inputdlg(promptCond, 'Order conditions');
            app.options.ConditionOrder = condOrd;
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        function app = networkActivityApp
            createComponents(app)
            s = settings;
            if ~hasGroup(s, 'networkActivity')
                % save the settings
                addGroup(s, 'networkActivity');
                addSetting(s.networkActivity, 'StillName', 'PersonalValue', 'DIC.');
                addSetting(s.networkActivity, 'Microscope', 'PersonalValue', 'Nikon A1');
                addSetting(s.networkActivity, 'Frequency', 'PersonalValue', 0);
                addSetting(s.networkActivity, 'RoiSize', 'PersonalValue', 5);
                addSetting(s.networkActivity, 'RoiShape', 'PersonalValue', 'Square');
                addSetting(s.networkActivity, 'PeakMinHeight', 'PersonalValue', 'MAD');
                addSetting(s.networkActivity, 'SigmaThr', 'PersonalValue', 2);
                addSetting(s.networkActivity, 'PeakMinProminance', 'PersonalValue', 0.1);
                addSetting(s.networkActivity, 'PeakMinDistance', 'PersonalValue', 1);
                addSetting(s.networkActivity, 'PeakMinDuration', 'PersonalValue', 2);
                addSetting(s.networkActivity, 'PeakMaxDuration', 'PersonalValue', 5);
                addSetting(s.networkActivity, 'DetectTrace', 'PersonalValue', 'Smooth');
                addSetting(s.networkActivity, 'LastPath', 'PersonalValue', pwd);
                addSetting(s.networkActivity, 'Registration', 'PersonalValue', false);
                addSetting(s.networkActivity, 'Reference', 'PersonalValue', 'Tyrode');
                addSetting(s.networkActivity, 'Detrending', 'PersonalValue', 'Moving median');
                addSetting(s.networkActivity, 'DetrendSize', 'PersonalValue', 100);
            end
            app.options.StillName = s.networkActivity.StillName.ActiveValue;
            app.options.Microscope = s.networkActivity.Microscope.ActiveValue;
            app.options.Frequency = s.networkActivity.Frequency.ActiveValue;
            app.options.RoiSize = s.networkActivity.RoiSize.ActiveValue;
            app.options.RoiShape = s.networkActivity.RoiShape.ActiveValue;
            app.options.PeakMinHeight = s.networkActivity.PeakMinHeight.ActiveValue;
            app.options.SigmaThr = s.networkActivity.SigmaThr.ActiveValue;
            app.options.PeakMinProminance = s.networkActivity.PeakMinProminance.ActiveValue;
            app.options.PeakMinDistance = s.networkActivity.PeakMinDistance.ActiveValue;
            app.options.PeakMinDuration = s.networkActivity.PeakMinDuration.ActiveValue;
            app.options.PeakMaxDuration = s.networkActivity.PeakMaxDuration.ActiveValue;
            app.options.DetectTrace = s.networkActivity.DetectTrace.ActiveValue;
            app.options.LastPath = s.networkActivity.LastPath.ActiveValue;
            app.options.Registration = s.networkActivity.Registration.ActiveValue;
            app.options.Reference = s.networkActivity.Reference.ActiveValue;
            app.options.Detrending = s.networkActivity.Detrending.ActiveValue;
            app.options.DetrendSize = s.networkActivity.DetrendSize.ActiveValue;
            if nargout == 0
                clear app
            end
        end
        
        function saveSettings(app)
            % Before closing save the settings
            s = settings;
            s.networkActivity.StillName.PersonalValue = app.options.StillName;
            s.networkActivity.Microscope.PersonalValue = app.options.Microscope;
            s.networkActivity.Frequency.PersonalValue = app.options.Frequency;
            s.networkActivity.RoiSize.PersonalValue = app.options.RoiSize;
            s.networkActivity.RoiShape.PersonalValue = app.options.RoiShape;
            s.networkActivity.PeakMinHeight.PersonalValue = app.options.PeakMinHeight;
            s.networkActivity.SigmaThr.PersonalValue = app.options.SigmaThr;
            s.networkActivity.PeakMinProminance.PersonalValue = app.options.PeakMinProminance;
            s.networkActivity.PeakMinDistance.PersonalValue = app.options.PeakMinDistance;
            s.networkActivity.PeakMaxDuration.PersonalValue = app.options.PeakMaxDuration;
            s.networkActivity.DetectTrace.PersonalValue = app.options.DetectTrace;
            s.networkActivity.LastPath.PersonalValue = app.options.LastPath;
            s.networkActivity.Registration.PersonalValue = app.options.Registration;
            s.networkActivity.Reference.PersonalValue = app.options.Reference;
        end
        
        function delete(app)
            % If the data is not saved, ask if it needs to be saved
            answer = questdlg('Do you want to save the data?', 'Save before closing');
            switch answer
                case 'Yes'
                    FileMenuSaveSelected(app);
                case 'No'
                    % Nothing to add
                case 'Cancel'
                    return
            end
            delete(app.Figure)
        end
    end
end