function [imSeq, images] = westheimer(stimset, video, timing)

    %% create base element
    % each matrix is of the same size
    % algorithm:
    
    video.minDim = min(video.width_pix, video.height_pix);
    
    stimset.baseElemSizePx = round(video.minDim/stimset.nBaseElems);
    stimset.pedestalSizePct = linspace(stimset.pedstalWidthPct(1), stimset.pedstalWidthPct(2), timing.nCoreSteps);
    stimset.probeLum = linspace(stimset.probeLumPct(1), stimset.probeLumPct(2), timing.nCoreSteps);

    timing.nUniqueFrames = timing.nCoreSteps*timing.figFramesPerCycle/...
        timing.sweepUpdateFrames;
    timing.nFramesPerStep = timing.nUniqueFrames/timing.nCoreSteps;
    
    
    frames = zeros( video.minDim,  video.minDim, 1, timing.nUniqueFrames);
    for s = 1:timing.nCoreSteps     
        frames(:, :, 1, (s - 1)*timing.nFramesPerStep +1:s*timing.nFramesPerStep) = mkSweepFrames(stimset, timing, video, s);
    end
    
    %% add ciolor dim
    images = uint8(frames);
    
    %% make imageSequence
    nCyclesPerStep = timing.framesPerBin/(timing.nFramesPerStep*timing.sweepUpdateFrames);
    imSeq = [];
    for s = 1:timing.nCoreSteps
        stepFrames = zeros(timing.framesPerBin/nCyclesPerStep, 1);
        stepFrames(1:timing.sweepUpdateFrames:end) = (s - 1)*timing.nFramesPerStep + 1:s*timing.nFramesPerStep;
        imSeq = cat(1, imSeq, repmat(stepFrames, [nCyclesPerStep 1]));
    end
    imSeq = uint32(imSeq);
end

function frames = mkSweepFrames(stimset, timing, video, nStep)
    % get sweep profile
    switch stimset.modulationVar
        case 'pedestal' 
            % pedestal grows
            pedestalSizePx = round(stimset.baseElemSizePx*stimset.pedestalSizePct(nStep))*ones(1, timing.nFramesPerStep);
            %probe updates every flickerFramesPerCycle
            probeLum = ones(1, timing.nCoreSteps);
            probeLum(1:timing.flickerFramesPerCycle:end) = 1 + stimset.probeLumPct(end);
            probeLum = stimset.probeLumBase.*probeLum;
        case 'probe'
            probeLum = stimset.probeLumBase*getModProfile(stimset.probeLum(nStep), timing.nFramesPerStep, stimset.modulationProfile);
            pedestalSizePx = round(ones(1, timing.nFramesPerStep)*stimset.baseElemSizePx*stimset.pedestalWidthBase);
    end
    
    frames = zeros(video.minDim, video.minDim, timing.nFramesPerStep);
    for nf = 1:timing.nFramesPerStep
        baseElem = makeSquareElement(stimset.baseElemSizePx, stimset.baseLum, pedestalSizePx(nf), stimset.pedestalLum, stimset.probeWidthPix, probeLum(nf));
        frames(:, :, nf) = repmat(baseElem, stimset.nBaseElems);
    end

end
function out = getModProfile(sweepValue, nPoints, modProfile)
    switch modProfile
        case 'sawtooth-on'
            out = linspace(0, sweepValue, nPoints);
        case 'sawtooth-off'
            out = linspace(sweepValue, 0, nPoints);
    end
end


