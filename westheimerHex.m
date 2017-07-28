    function [imSeq, images] = westheimerHex(stimset, videoMode, timing)
        
        %% create base element
        % each matrix is of the same size
        % algorithm:
        
        videoMode.minDim = min(videoMode.width_pix, videoMode.height_pix);
        videoMode.width_deg = 2 * atand( (videoMode.width_cm/2)/videoMode.viewDistCm );
        videoMode.height_deg = 2 * atand( (videoMode.height_cm/2)/videoMode.viewDistCm );
        videoMode.pix2arcmin = ( videoMode.width_deg * 60 ) / videoMode.width_pix;
        
        timing.nFramesPerStep = nUniqueFramesPerStep(stimset, timing);
        
        timing.UniqueSteps = timing.nCoreSteps;
        if (strcmp(stimset.sweepType, 'Fixed'))
            timing.UniqueSteps = 1;
        end
        
        timing.nUniqueFrames = timing.UniqueSteps*timing.nFramesPerStep;

        stimset.maxLuminance = 255;
        % hexagons
        stimset.nSides = 6;
        
        stimset.baseSizePix = floor(stimset.baseSizeAmin/videoMode.pix2arcmin);
        stimset.pedestalSizePix = floor(stimset.pedestalSizePct*stimset.baseSizePix/100);
        stimset.probeSizePix = floor(stimset.probeSizePct*stimset.baseSizePix/100);
		
		stimset.pedestalLuminance = floor(stimset.maxLuminance*stimset.pedestalLum/100);
        stimset.backgroundLuminance = floor(stimset.maxLuminance*stimset.bgrLum/100);
        stimset.probeLuminanceMax = floor(stimset.pedestalLuminance*(1 + stimset.probeLuminance/100));
        
        frames = mkTrialSequence(stimset, videoMode, timing);
        
       
        %% add color dim
        images = uint8(frames);
        
        %% make imageSequence
        imSeq = [];
        updateEveryFrame = timing.nFramesPerStep;
        totalRepeats = timing.framesPerStep/updateEveryFrame;
        uniqueSteps = linspace(1, timing.UniqueSteps, timing.nCoreSteps);
        
        for s = 1:timing.nCoreSteps
            stepFrames = (1 + timing.nFramesPerStep*(uniqueSteps(s) - 1):uniqueSteps(s)*timing.nFramesPerStep)';
            imSeq = cat(1, imSeq, repmat(stepFrames, [totalRepeats 1]));
        end
		preludeSeq = [];
		postludeSeq = [];
		if (timing.nPreludeBins)
			% repeat first step and last step
			preludeSeq = imSeq(1:timing.nPreludeFrames);
			postludeSeq = imSeq(end - timing.nPreludeFrames + 1:end);
		end
		
        imSeq = [preludeSeq; uint32(imSeq); postludeSeq];
    end
    
    function frames = mkTrialSequence(stimset, video, timing)
    
        
        %max number of hexagons after magnification effect:    
        nMaxHex = getMaxHex(video.minDim, stimset.mfactor, stimset.baseSizePix, video.pix2arcmin/60);
    
        %start at top left corner and map lower quarter (IV Quadrant))         
        x0 = 0;
        y0 = 0;
        % hex centers
        hexCenters = getHexCenters(x0, y0, stimset.baseSizePix, nMaxHex);
        quarter = zeros(.5*video.minDim); 
        pedestalCenters = cell2mat(hexCenters);
        
        %% remove hexagons lying on the X OR Y axis       
        outIdxPd = (pedestalCenters(:, 1).*pedestalCenters(:, 2))>0;
        pedestalCenters = pedestalCenters(outIdxPd, :);
        pedestalCenters = [[0, 0]; pedestalCenters];
        %% Probe map (same for the entire sequence)
        
        probeActivationTable = probeLookupTable(size(hexCenters, 1), stimset.probeTable);
        probeCenters = cell2mat(hexCenters(probeActivationTable >0));

        outIdxPr = (probeCenters(:, 1).*probeCenters(:, 2))>0;        
        probeCenters = probeCenters(outIdxPr, :);

        probeHexMap_q4 = makeMask_CM(probeCenters, stimset.probeSizePix, ...
            stimset.nSides, size(quarter), stimset.mfactor, video.pix2arcmin/60);            
        
       
        hexProbeMap = mkWhole(probeHexMap_q4);
        
        sfMask = applySubfieldMask(hexProbeMap, stimset.UpDown, stimset.LeftRight); 
        hexProbeMap = hexProbeMap.*sfMask;
        
        %% 3. Constuct frames
        
        sweepModProfile = getModProfile(0, 1, timing.nFramesPerStep, stimset.modType);
        baseMap = stimset.backgroundLuminance*ones(size(hexProbeMap));
        
        hexProbeFrames = repmat(hexProbeMap, [1 1 timing.nFramesPerStep]);

        hexFramesCell = cell(timing.UniqueSteps, 1);
        for s = 1:timing.UniqueSteps
            pedestalHexMap_q4 = makeMask_CM(pedestalCenters, stimset.pedestalSizePix(s), ...
                stimset.nSides, size(quarter), stimset.mfactor, video.pix2arcmin/60);    
            hexPedestalMapStep = mkWhole(pedestalHexMap_q4);
            
            
            pedestalFramesStep = (stimset.pedestalLuminance - stimset.backgroundLuminance)...
                *repmat(hexPedestalMapStep, [1 1 timing.nFramesPerStep]);
            
            probeLuminanceStep = 2.5*sweepModProfile.*(stimset.probeLuminanceMax(s) - stimset.pedestalLuminance);
            
            probeFramesStep = bsxfun(@times, hexProbeFrames,...
                reshape(probeLuminanceStep, 1, 1, numel(probeLuminanceStep)));
            
            baseFramesStep = repmat(baseMap, [1 1 timing.nFramesPerStep]);
            
            hexFramesCell{s} = uint8(pedestalFramesStep + baseFramesStep + probeFramesStep);
        end
        
        hexFrames =  cell2mat(hexFramesCell);
            
%         % Place fixation cross        
%         baseCrossFrame = .5*ones(size(hexFrames, 1), size(hexFrames, 2));
% 		offset = [0, 0];
% 		height = 30;
%         crossFrame = addCross(baseCrossFrame, height, offset);
% 		crossFrame(crossFrame <= 0.5) = 0;
% 		crossFrames = uint8(255*repmat(crossFrame, [1 1  size(hexFrames, 3)]));
%         frames(:, :, 1, :) = crossFrames + hexFrames;
		frames(:, :, 1, :) = hexFrames;
    end
    function nHex = getMaxHex(height, k, r, pix2deg)
        % max visual field degree, half
        max_Xdeg = (0.5*height)*pix2deg;
    
        %what point will project to max visual field degree, half
        max_Xdeg_new = 0.5*(-k + sqrt(k^2 + 4*k*max_Xdeg))/pix2deg;
        nHex = round(max_Xdeg_new/(2*r));
    end

    function xyPix_new = getMaxXY(xyCoord, k, r, pix2deg)
        % max visual field degree, half
        xyDeg0 = xyCoord*pix2deg;
        xyDeg1 = (xyCoord + r)*pix2deg;

        %what point will project to max visual field degree, half
        xy_Deg0_new = xyDeg0*(1 + xyDeg0/k);                
        xy_Deg1_new = xyDeg1*(1 + xyDeg1/k);               
        
        xyPix_new = (xy_Deg1_new - xy_Deg0_new)/pix2deg;
    end
    
    function hexCenter = getHexCenters(x0, y0, rp, nMaxHex)
        
        % large radius
        %rp = round(Rp*sqrt(3)/2);
        
        Rp = round(rp*2/sqrt(3));        
        hexCenter = cell(nMaxHex + 1, 1);
        hexCenter{1} = [x0 y0];  
        for h = 1:nMaxHex;       
            %diagonal
            x_d = x0 + rp*(2*h:-1:h)';
            y_d = y0 + Rp*(0:1.5:1.5*h)';
        
            %horizontal
            xh_num = floor(h/2);
            xh_start = rem(h, 2);
            xh_end = xh_start + 2*(xh_num - 1);
            x_h = x0 + rp*(xh_start:2:xh_end)';
            y_h = y0 + Rp*ones(numel(x_h), 1)*3/2*h;
            hexCenter{h + 1} = [[x_d y_d]; [x_h y_h]];
        end
    end

    
    function [XM, YM] = someShapeMagnified(center, radius, n, k, pix2deg)
    
        t = linspace(0, 2*pi, n);
        r = radius*ones(1, n);
    
        X = r.*sin(t) + center(1);
        Y = r.*cos(t) + center(2);
    
        xD = X*pix2deg;
        yD = Y*pix2deg;
    
        [theta, rho] = cart2pol(xD, yD);
            
        s2 = rho.*(1 + rho/k);

        [XM_deg, YM_deg] =  pol2cart(theta, s2);
    
        XM = XM_deg/pix2deg;
        YM = YM_deg/pix2deg;
    end
%     function [X, Y] = someShape(center, radius, n, k, pix2deg)
%     
%         t = linspace(0, 2*pi, n);
%         r = radius*ones(1, n);
%     
%         X = r.*sin(t) + center(1);
%         Y = r.*cos(t) + center(2);    
%     end

    function maskOut = makeMask_CM(center, R, n_sides, msize, mfactor, pix2deg)
        maskOut = zeros(msize);
        % center coordiantes kth row = kth polyshape, col1 = x, col2 = y
        for c = 1:size(center, 1)
            [xr, yr] = someShapeMagnified(center(c, :), R, n_sides + 1, mfactor, pix2deg);
            maskOut = maskOut + poly2mask(floor(xr), floor(yr), msize(1), msize(2));
        end 
    end
    function maskOut = makeMask(center, R, n_sides, msize, mfactor, pix2deg)
        maskOut = zeros(msize);
        % center coordiantes kth row = kth polyshape, col1 = x, col2 = y
        for c = 1:size(center, 1)
            [xr, yr] = someShape(center(c, :), R, n_sides + 1, mfactor, pix2deg);
            maskOut = maskOut + poly2mask(floor(xr), floor(yr), msize(1), msize(2));
        end 
    end
    
    %% Sector map 
    function probeTable = probeLookupTable(nCenters, sector)
        
        nSectors = 3;
        % remove the center element
        sectorWidth = round((nCenters - 1)/nSectors);
        sectorOverlap = 1;
        probeTable = zeros(nCenters, 1);
        switch sector
            case 'f1'                
                fStart = 1;
                fEnd = sectorWidth + sectorOverlap + 1;
            case 'p1'
                % extra 1 for the center probe (first in the map)
                fStart = 1 + sectorWidth + sectorOverlap;
                fEnd = fStart + sectorWidth;
            case 'p2'
                 % extra 1 for the center probe (first in the map)               
                fStart = 1 + 1 + 2*sectorWidth - sectorOverlap;
                fEnd = nCenters;
                
            case 'All'
                fStart = 1;
                fEnd = nCenters;
        end
        probeTable(fStart:fEnd) = 1;
    end
    
    
    function out = mkWhole(quarter) 
        hexMap_tophalf= [quarter(end:-1:1, end:-1:1) quarter(end:-1:1, :)];
        out = [hexMap_tophalf; hexMap_tophalf(end:-1:1, :)];                
    end
    
    function out = applySubfieldMask(frame, InfSup, TempNas)
    
        quarterMask = ones(0.5*size(frame));
        is = zeros(2, 2);
        tn = zeros(2, 2);
        switch InfSup
            case 'Upper'
                is(1, :) = 1;
            case 'Lower'
                is(2, :) = 1;
            case 'All'
                is(:, :) = 1;
        end
    
        switch TempNas
            case 'Left'
                tn(:, 1) = 1;
            case 'Right'
                tn(:, 2) = 1;
            case 'All'
                tn(:,:) = 1;  
        end
        
        q = is.*tn;
        out = [q(1, 1)*quarterMask, q(1, 2)*quarterMask; q(2, 1)*quarterMask, q(2, 2)*quarterMask];
    end
    
    %% Sweep frames   
    
    function nf = nUniqueFramesPerStep(stimset, timing)
        switch(stimset.modType)
            case 'None'
                nf = timing.updateFramesPerCycle;
            case 'Square'
                nf = timing.updateFramesPerCycle;
            case {'Sawtooth-on', 'Sawtooth-off'}
                nf = timing.updateFramesPerCycle;
        end
        
    end
    function out = getModProfile(firstVal, lastVal,  nPoints, modProfile)
        switch modProfile
            case 'None'
                out = firstVal*ones(nPoints, 1);
            case 'Square'
                out = [firstVal*ones(round(.5*nPoints), 1); lastVal*ones(nPoints - round(.5*nPoints), 1)];
            case 'Sawtooth-on'
                out = linspace(firstVal, lastVal, nPoints);
            case 'Sawtooth-off'
                out = linspace(lastVal, firstVal, nPoints);                
        end
    end

	function out = addCross(image, height, offset)
		
		%default value of cross liner;
		outer = 0.5;
		
		if (size(image, 3) > 1)
			outer = mean(mean(rgb2gray(image)));
		else
			outer = mean(image(:));
		end
		% zero offset = > place cross exactly on image center
		if (nargin < 3)
			offset = [0, 0];
		end;
		
		% outer parameters: mean lum, 1 pix off center
		outer_thick = 1;
		outer_code = 1;
		% inner parameters: white, 1 pix off center
		inner_code = 0.5;
		inner_thick = 1;
		
		
		dummy = zeros(size(image, 1), size(image, 2));
		center0 = floor([size(image, 1), size(image, 2)]*0.5);
		
		center = center0 + offset;
		%draw inner
		dummy(center(1) - inner_thick:center(1) + inner_thick, center(2) - height:center(2) + height) = inner_code;
		dummy(center(1) - height:center(1) + height, center(2) - inner_thick:center(2) + inner_thick) = inner_code;
		
		%outer boundaries
		outer_boundaryV = [center(1) - height - 1, center(1) - inner_thick - 1, center(1) + inner_thick + 1, center(1) + height + 1];
		outer_boundaryH = [center(2) - height - 1, center(2) - inner_thick - 1, center(2) + inner_thick + 1, center(2) + height + 1];
		
		%draw vertical outer
		dummy(outer_boundaryV(1):outer_boundaryV(2), outer_boundaryH(2) - 2*outer_thick:outer_boundaryH(2)) = outer_code;
		dummy(outer_boundaryV(1):outer_boundaryV(2), outer_boundaryH(3):outer_boundaryH(3) + 2*outer_thick) = outer_code;
		
		dummy(outer_boundaryV(2):outer_boundaryV(3), outer_boundaryH(1) -  2*outer_thick:outer_boundaryH(1)) = outer_code;
		dummy(outer_boundaryV(2):outer_boundaryV(3), outer_boundaryH(4):outer_boundaryH(4) + 2*outer_thick) = outer_code;
		
		dummy(outer_boundaryV(3):outer_boundaryV(4), outer_boundaryH(2)- 2*outer_thick:outer_boundaryH(2)) = outer_code;
		dummy(outer_boundaryV(3):outer_boundaryV(4), outer_boundaryH(3):outer_boundaryH(3) + 2*outer_thick) = outer_code;
		
		%draw horizontal outer
		dummy(outer_boundaryV(1) - 2*outer_thick:outer_boundaryV(1), outer_boundaryH(2):outer_boundaryH(3)) = outer_code;
		dummy(outer_boundaryV(4):outer_boundaryV(4) + 2*outer_thick, outer_boundaryH(2):outer_boundaryH(3)) = outer_code;
		
		dummy(outer_boundaryV(2) - 2*outer_thick:outer_boundaryV(2), outer_boundaryH(1):outer_boundaryH(2)) = outer_code;
		dummy(outer_boundaryV(3):outer_boundaryV(3) + 2*outer_thick, outer_boundaryH(1):outer_boundaryH(2)) = outer_code;
		
		dummy(outer_boundaryV(2) - 2*outer_thick:outer_boundaryV(2), outer_boundaryH(3):outer_boundaryH(4)) = outer_code;
		dummy(outer_boundaryV(3):outer_boundaryV(3) + 2*outer_thick, outer_boundaryH(3):outer_boundaryH(4)) = outer_code;
		
		dummy = repmat(dummy,[1,1,size(image,3)]);
		out = image;
		out(dummy == 0.5) = 0.6;
		out(dummy == 1) = outer;
    end
