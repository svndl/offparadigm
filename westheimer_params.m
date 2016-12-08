function [stimset, video, timing] = westheimer_params

    %% video settings
    
    video.extent = 'square';
    video.width_pix = 2560;
    video.height_pix = 1600;
    
    %% stimset settings
    
    stimset.baseLum = 0;
    stimset.nBaseElems = 10;
    stimset.probeWidthPix = 10;
    
    %percent of base size
    stimset.pedestalWidthBase = .07;
    stimset.pedstalWidthPct = [.07, .2];
    stimset.pedestalLum = 127;
    stimset.modulationProfile = 'sawtooth-on';
    
    % swept, use an array
    stimset.modulationVar = 'pedestal';
    stimset.probeLumBase = stimset.pedestalLum;
    stimset.probeLumPct = [0, .2];
    
    
    %% timing
    timing.nCoreSteps = 10;
    timing.framesPerBin = 60;
    timing.sqrFramesPerCycle = 6;
    timing.figFramesPerCycle = 60;
    timing.nCoreFrames = 600;
end

