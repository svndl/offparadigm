function pmf_RandomDotsStereoMotion( varargin )
%
%   Shows a correlated figure region which potentially modulates in depth
%   over a background whose correlation level can be any value from 0-1.
%
%	xDiva Matlab Function paradigm

%	"pmf_" prefix is not strictly necessary, but helps to identify
%	functions that are intended for this purpose.

%   Each Matlab Function paradigm will have its own version of this file

%	First argument is string for selecting which subfunction to use
%	Additional arguments as needed for each subfunction

%	Each subfunction must conclude with call to "assignin( 'base', 'output', ... )",
%	where value assigned to "output" is a variable or cell array containing variables
%	that xDiva needs to complete the desired task.

if nargin > 0, aSubFxnName = varargin{1}; else error( 'pmf_RandomDotsStereoMotion.m called with no arguments' ); end

% these next three are shared by nested functions below, so we create
% them in this outermost enclosing scope.
definitions = MakeDefinitions;
parameters = {};
timing = {};
videoMode = {};

% some useful functional closures...
CFOUF = @(varargin) cellfun( varargin{:}, 'uniformoutput', false );
AFOUF = @(varargin) arrayfun( varargin{:}, 'uniformoutput', false );

%lambda functions
PVal = @( iPart, x ) ParamValue( num2str(iPart), x );
PVal_S = @(x) ParamValue( 'S', x );



% PVal_B = @(x) ParamValue( 'B', x );
% PVal_1 = @(x) ParamValue( 1, x );
% PVal_2 = @(x) ParamValue( 2, x );

aaFactor = 8;
maxDispAmin = 70;
ppath = setPathParadigm;
screenRedGunData = load(fullfile(ppath.info, 'SonyTV_RedGunValues.txt'));

try
    switch aSubFxnName
        case 'GetDefinitions', GetDefinitions;
        case 'ValidateParameters', ValidateParameters;
        case 'MakeMovie', MakeMovie;
    end
catch tME
    errLog = fopen(fullfile(ppath.log, 'ErrorLog.txt'), 'a+');
    display(tME.message);
    for e = 1: numel(tME.stack)
        fprintf(errLog, ' %s ', tME.stack(e).file);
        fprintf(errLog, ' %s ', tME.stack(e).name);
        fprintf(errLog, ' %d\n', tME.stack(e).line);
    end
    fclose(errLog);
    rethrow( tME ); % this will be caught by xDiva for runtime alert message
end

    function rV = ParamValue( aPartName, aParamName )
        % Get values for part,param name strings; e.g "myViewDist = ParamValue( 'S', 'View Dist (cm)' );"
        tPart = parameters{ ismember( { 'S' 'B' '1' '2' }, {aPartName} ) }; % determine susbscript to get {"Standard" "Base" "Part1" "Part2"} part cell from parameters
        rV = tPart{ ismember( tPart(:,1), { aParamName } ), 2 }; % from this part, find the row corresponding to aParamName, and get value from 2nd column
    end
    function rV = GetParamArray( aPartName, aParamName )
        
        % For the given part and parameter name, return an array of values
        % corresponding to the steps in a sweep.  If the requested param is
        % not swept, the array will contain all the same values.
        
        % tSpatFreqSweepValues = GetParamArray( '1', 'Spat Freq (cpd)' );
        
        % Here's an example of sweep type specs...
        %
        % definitions{end-2} =
        % 	{
        % 		'Fixed'         'constant'   { }
        % 		'Contrast'      'increasing' { { '1' 'Contrast (pct)' } { '2' 'Contrast (pct)' } }
        % 		'Spat Freq'      'increasing' { { '1' 'Spat Freq (cpd)' } { '2' 'Spat Freq (cpd)' } }
        % 	}
        
        T_Val = @(x) timing{ ismember( timing(:,1), {x} ), 2 }; % get the value of timing parameter "x"
        tNCStps = T_Val('nmbCoreSteps');
        tSweepType = PVal_S('Sweep Type');
        
        % we need to construct a swept array if any of the {name,value} in definitions{5}{:,3}
        
        [ ~, tSS ] = ismember( tSweepType, definitions{end-2}(:,1) ); % the row subscript in definitions{5} corresponding to requested sweep type
        % determine if any definitions{5}{ tSS, { {part,param}... } } match arguments tPartName, tParamName
        IsPartAndParamMatch = @(x) all( ismember( { aPartName, aParamName }, x ) );
        tIsSwept = any( cellfun( IsPartAndParamMatch, definitions{end-2}{tSS,3} ) ); % will be false for "'Fixed' 'constant' { }"
        
        if ~tIsSwept
            rV = ones( tNCStps, 1 ) * ParamValue(  aPartName, aParamName );
        else
            tStepType = PVal_S('Step Type');
            tIsStepLin = strcmpi( tStepType, 'Lin Stair' );
            tSweepStart = PVal_S('Sweep Start');
            tSweepEnd = PVal_S('Sweep End');
            if tIsStepLin
                rV = linspace( tSweepStart, tSweepEnd, tNCStps )';
            else
                rV = logspace( log10(tSweepStart), log10(tSweepEnd), tNCStps )';
            end
        end
        
    end

    function rV = MakeDefinitions
        % for "ValidateDefinition"
        % - currently implementing 'integer', 'double', 'nominal'
        % - types of the cells in each parameter row
        % - only "standard" type names can be used in the "type" fields
        % - 'nominal' params should have
        %       (a) at least one item
        %		(b) value within the size of the array
        % - all other params (so far) should have empty arrays of items
        
        rV = { ...
            
        % - Parameters in part_S must use standard parameter names
        % - 'Sweep Type' : at least 'Fixed' sweep type must be defined as first item in list
        % - 'Modulation' : at least 'None' modulation type must be defined
        % - 'Step Type'  : at least 'Lin Stair' type must be defined,
        %                  first 4 step types are reserved, custom step types can only be added after them
        
        % "Standard" part parameters - common to all paradigms, do not modify names.
        {
        'View Dist (cm)'	100.0	        'double'  {}
        'Mean Lum (cd)'     1.0	            'double'  {} % under default calibration, this is (0,0,0)
        'Fix Point'         'None'	        'nominal' { 'None' 'Cross' }
        'Sweep Type'        'Pedestal Size' 'nominal' { 'Fixed' 'Pedestal Size' 'Probe Contrast' }
        'Step Type'         'Lin Stair'	    'nominal' { 'Lin Stair' 'Log Stair' }
        'Sweep Start'       0.0	            'double'  {}
        'Sweep End'         1.0             'double'  {}
        'Modulation'        'None'          'nominal' { 'None' 'Sawtooth' 'Square' } %   
        }
        %
        
        % "Base" part parameters - paradigm specific parameters that apply to unmodulated parts of the stimulus
        {
        'ModInfo'                   0.0       'integer' {}
        'Squares per row'	        10.0      'integer'	{}
        'Bgr lum (% max)'           10        'integer' {} 
        'Pedestal lum (% max)'      50        'integer' {} 
        'Pedestal Size (amin)'      5.0       'double'  {}        
        }
        
        % "Part1" - parameters that apply to part of stimulus that carries first frequency tag.
        % "Cycle Frames" must be first parameter
        {
        'Cycle Frames'           2.0         'integer'	{}  % framerate(Hz)/stimFreq(Hz) Modulation frequency every other frame (2 Hz)
        'Contrast (pct)'	     20.0        'double'	{}
        'Probe Polarity'         'positive'  'nominal'  {'positive', 'negative'}
        'Probe size (amin)'      5.0         'double'   {}         
        }
        
        
        % Sweepable parameters
        % The cell array must contain as many rows as there are supported Sweep Types
        % 1st column (Sweep Types) contains Sweep Type as string
        % 2nd column (Stimulus Visiblity) contains one of the following strings,
        % indicating how stimulus visibility changes when corresponding swept parameter value increases:
        %   'constant' - stimulus visibility stays constant
        %   'increasing' - stimulus visibility increases
        %   'decreasing' - stimulus visibility decreases
        % 3rd column contains a single-row cell array of pairs, where each pair is a single-row cell
        % array of 2 strings: { Part name, Parameter name }
        
        % If sweep affects only one part, then you only need one
        % {part,param} pair; if it affects both parts, then you need both
        % pairs, e.g. for "Contrast" and "Spat Freq" below
        
        {
        'Fixed'			   'constant'   { }
        'Pedestal Size'    'constant' { { 'B' 'Pedestal Size (amin)' } }
        'Probe Contrast'   'increasing' { { '1' 'Contrast (pct)' } }
        }        
        
        % ModInfo information
        % The cell array must contain as many rows as there are supported Modulations
        % 1st column (Modulation) contains one of the supported Modulation typs as string
        % 2nd column contains the name of the ModInfo parameter as string
        % 3rd column (default value) contains default value of the ModInfo
        % parameter for this Modulation
        {        
         'None'             'ModInfo'	       0.0
         'Sawtooth'         'ModInfo'          0.0
         'Square'           'ModInfo'          0.0
         }
        % check direction
        % Required by xDiva, but not by Matlab Function
        {
        'Version'					1
        'Adjustable'				true
        'Needs Unique Stimuli'		false % ###HAMILTON for generating new stimuli every time
        'Supports Interleaving'		false
        'Part Name'                 { 'Probe' }
        'Frame Rate Divisor'		{ 2 } % {even # frames/cycle only, allows for odd-- makes sense for dot update}
        'Max Cycle Frames'			{ 60 } % i.e. -> 0.5 Hz, 10 Hz
        'Allow Static Part'			{ true }
        }
        };
    end

    function GetDefinitions
        assignin( 'base', 'output', MakeDefinitions );
    end

    function ValidateParameters
        % xDiva invokes Matlab Engine command:
        
        % pmf_<subParadigmName>( 'ValidateParameters', parameters, timing, videoMode );
        % "parameters" here is an input argument. Its cellarray hass the
        % same structure as "defaultParameters" but each parameter row has only first two
        % elements
        
        % The "timing" and "videoMode" cellarrays have the same row
        % structure with each row having a "name" and "value" elements.
        
        
        [ parameters, timing, videoMode ] = deal( varargin{2:4} );
 
        %% get video system info
                
        VMVal = @(x) videoMode{ ismember( videoMode(:,1), {x} ), 2 };
        
        width_pix = VMVal('widthPix');
        height_pix = VMVal('heightPix');
        
        width_cm = VMVal('imageWidthCm');
        viewDistCm = PVal('S','View Dist (cm)');
        width_deg = 2 * atand( (width_cm/2)/viewDistCm );
        pix2arcmin = ( width_deg * 60 ) / width_pix;
        
        probeSizeAmin = PVal('1', 'Probe size (amin)');
        probeSizePix = probeSizeAmin/pix2arcmin;
        pedestalSizeAmin = PVal('B', 'Pedestal Size (amin)');
        pedestalSizePix = pedestalSizeAmin/pix2arcmin;        
        nSquares = PVal('B','Squares per row');
        squareSizePix = round(width_pix/nSquares);
          
        validationMessages = {};

        
        % size validation
        ValidateProbeSize;
        isSwept = ~strcmp(PVal('S', 'Sweep Type'), 'Fixed');
        if (isSwept)
            % sweep range validation
            ValidateSweepRange;
        else
            ValidatePedestalSize;
            ValidateProbeLum;
        end
        ValidateBackgroundSize;
        %luminance validation
        ValidateBackgroundLum;
        ValidatePedestalLum;
        
        
        % Standard routine
        parametersValid = isempty( validationMessages );
        output = { parametersValid, parameters, validationMessages };
        assignin( 'base', 'output', output );

        
        % min and max size (1 pix - background square region) 
        function ValidateProbeSize
            corrected = 0;
            if mod(probeSizePix, 1) %decimal
                probeSizePix = round(probeSizePix);
                corrected = 1;
                msg = 'not integer';
            end
            if probeSizePix < 1
                probeSizePix = 1;
                corrected = 1;
                msg = 'too small';
            end
            if (probeSizePix > squareSizePix)
                probeSizePix = .8*squareSizePix;
                corrected = 1;
                msg = 'too large';
            end
            if (corrected)
                 probeSizeAmin = probeSizePix*pix2arcmin;
                CorrectParam('1', 'Probe size (amin)', probeSizeAmin);
                AppendVMs(sprintf(...
                    'Probe size is %s corrected to nearest possible value: %3.4f amin.',...
                    msg, probeSizeAmin));
            end
               
        end
        
        % min and max pedestal (1 pix - background square region )
        function ValidatePedestalSize
            corrected = 0;            
            if mod(pedestalSizePix, 1) %decimal
                pedestalSizePix = round(pedestalSizePix);                
                corrected = 1;
                msg = 'not integer';
            end
            if pedestalSizePix < 1
                pedestalSizePix = 1;
                msg = 'too small';
                corrected = 1;                
            end
            if (pedestalSizePix > squareSizePix)
                pedestalSizePix = .9*squareSizePix;
                msg = 'too large';
                corrected = 1;                
            end
            if (corrected)
                pedestalSizeAmin = pedestalSizePix*pix2arcmin;
                CorrectParam('B', 'Pedestal Size (amin)', pedestalSizeAmin);
                AppendVMs(sprintf(...
                    'Pedestal size %s, corrected to nearest possible integer pixel value: %3.4f amin.',...
                    msg, pedestalSizeAmin));
            end
                
        end
        
        % limit the number of square elements 
        function ValidateBackgroundSize            
            %square size should be at least 2 pixel + probe + base            
            minSquareSize = probeSizePix + pedestalSizePix + 5;
            
            if (squareSizePix < minSquareSize)
                nSquares_new = floor(width_pix/minSquareSize);
                CorrectParam('B','Squares per row', nSquares_new);
                AppendVMs(sprintf(...
                    'Square size is too small, corrected to nearest possible value: %d ',...
                    nSquares_new));
            end
                
        end
        % max luminance is pedestal lum*( 1 + % ) 
        function ValidateProbeLum
        end

        % ValidateBackgroundLum?
        function ValidateBackgroundLum()
        end
        % ValidatePedestalLum?
        function ValidatePedestalLum()
        end
        
        % valid sweep range 
        function ValidateSweepRange()
            sweepType = PVal('S','Sweep Type');
            minSweep = PVal('S', 'Sweep Start');
            maxSweep = PVal('S', 'Sweep End');
                                    
            switch (sweepType)
                case 'Pedestal Size'
                    minPedestalSizePix = minSweep/pix2arcmin;
                    maxPedestalSizePix = maxSweep/pix2arcmin;
                    if (minPedestalSizePix < 1)
                        minPedestalSizePix = 1;
                        minSweep = minPedestalSizePix*pix2arcmin;
                        CorrectParam('S','Sweep Start', minSweep);
                        AppendVMs(sprintf( ...
                            'Sweep start is too small, corrected to nearest possible value: %3.4f', minSweep));
                    end
                    if (2*maxPedestalSizePix > squareSizePix)
                        maxPedestalSizePix = round(.9*squareSizePix/2);
                        
                        maxSweep = maxPedestalSizePix*pix2arcmin;
                        CorrectParam('S','Sweep End', maxSweep);
                        AppendVMs(sprintf( ...
                            'Sweep start is too small, corrected to nearest possible value: %3.4f', maxSweep));
                    end
                case 'Probe Luminance'

                    
            end
      
        end
                
        %% system function
        
        function CorrectParam( aPart, aParam, aVal )
            tPartLSS = ismember( { 'S' 'B' '1' '2' }, {aPart} );
            tParamLSS = ismember( parameters{ tPartLSS }(:,1), {aParam} );
            parameters{ tPartLSS }{ tParamLSS, 2 } = aVal;
        end
        
        function AppendVMs(aStr), validationMessages = cat(1,validationMessages,{aStr}); end       
    end

    function MakeMovie
        % ---- GRAB & SET PARAMETERS ----
        [ parameters, timing, videoMode, trialNumber ] = deal( varargin{2:5} );
        save('pmf_Westheimer_MakeMovie.mat', 'parameters', 'timing', 'videoMode');
        TRVal = @(x) timing{ ismember( timing(:,1), {x} ), 2 };
        VMVal = @(x) videoMode{ ismember( videoMode(:,1), {x} ), 2 };
        
        needsUnique = definitions{end}{3,2};
        needsImFiles = true;
        preludeType = {'Dynamic', 'Blank', 'Static'};
        % timing/trial control vars
        stimsetTiming.nCoreSteps = TRVal('nmbCoreSteps');
        stimsetTiming.nCoreBins = TRVal('nmbCoreBins');
        stimsetTiming.nPreludeBins = TRVal('nmbPreludeBins');
        stimsetTiming.framesPerStep = TRVal('nmbFramesPerStep');
        stimsetTiming.framesPerBin = TRVal('nmbFramesPerBin');
        stimsetTiming.preludeType = preludeType{1 + TRVal('preludeType')};
        stimsetTiming.isBlankPrelude = stimsetTiming.preludeType == 1;
        stimsetTiming.nCoreFrames = stimsetTiming.framesPerStep * stimsetTiming.nCoreSteps;
        stimsetTiming.nPreludeFrames = stimsetTiming.nPreludeBins * stimsetTiming.framesPerBin;
        stimsetTiming.nTotalFrames = 2 * stimsetTiming.nPreludeFrames + stimsetTiming.nCoreFrames;
        stimsetTiming.updateFramesPerCycle = PVal(1,'Cycle Frames'); % part 1 = Figure
        
        % screen vars
        video.width_pix = VMVal('widthPix');
        video.height_pix = VMVal('heightPix');
        video.width_cm = VMVal('imageWidthCm');
        video.height_cm = VMVal('imageHeightCm');
        video.frameRate = VMVal('nominalFrameRateHz');
        video.minLuminanceCd = VMVal('minLuminanceCd');
        video.maxLuminanceCd = VMVal('maxLuminanceCd');
        video.meanLuminanceCd = VMVal('meanLuminanceCd');
        video.meanLuminanceBitmap = VMVal('meanLuminanceBitmapValue');
        video.gammaTableCapacity = VMVal('gammaTableCapacity');
        
        video.viewDistCm = PVal('S','View Dist (cm)');
        video.aaFactor = aaFactor;
                
        % stim vars
        stimset.sweepType = PVal('S','Sweep Type');
        stimset.isSwept = ~strcmp(stimset.sweepType,'Fixed');  
        stimset.modType = PVal('S', 'Modulation');
        stimset.nBaseElems = PVal('B', 'Squares per row');
        stimset.bgrLum = PVal('B', 'Bgr lum (% max)' );
        
        stimset.pedestalLum = PVal('B', 'Pedestal lum (% max)');
        stimset.pedestalSizeAmin = GetParamArray('B', 'Pedestal Size (amin)');
        
        stimset.probeSizeAmin = PVal('1', 'Probe size (amin)');
        stimset.probeLuminance = GetParamArray('1', 'Contrast (pct)');
        stimset.probePolarity =  GetParamArray('1', 'Probe Polarity');      
        
        save('WestheimerInput.mat', 'stimset', 'video', 'stimsetTiming');
        
        [rImSeq, rIms] = westheimer(stimset, video, stimsetTiming);
        save('WestheimerOutput.mat', 'rImSeq', 'rIms');
                
        isSuccess = true;
        output = { isSuccess, rIms, cast( rImSeq, 'int32') }; % "Images" (single) and "Image Sequence" (Int32)
        %clear rIms
        assignin( 'base', 'output', output )
    end
    %%%%%%%%% STIMULUS GENERATION PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [imSeq, images] = westheimer(stimset, video, timing)
        
        %% create base element
        % each matrix is of the same size
        % algorithm:
        
        video.minDim = min(video.width_pix, video.height_pix);
        video.width_deg = 2 * atand( (video.width_cm/2)/video.viewDistCm );
        video.height_deg = 2 * atand( (video.height_cm/2)/video.viewDistCm );
        video.pix2arcmin = ( video.width_deg * 60 ) / video.width_pix;
        
        timing.nFramesPerStep = nUniqueFramesPerStep(stimset, timing);
        timing.nUniqueFrames = timing.nCoreSteps*timing.nFramesPerStep;
        
        %% recalculate base parameters (size, luminance)
        %     switch stimset.probePolarity
        %         case 'positive'
        %             stimset.probeContrastSign = 1;
        %         case 'negative'
        %             stimset.probeContrastSign = -1;
        %     end
        %stimset.maxLuminance = round(video.maxLuminanceCd*video.meanLuminanceBitmap/video.meanLuminanceCd);
        stimset.maxLuminance = 255;
        stimset.baseElemSizePx = round(video.minDim/stimset.nBaseElems);
        stimset.probeSizePix = round(stimset.probeSizeAmin/video.pix2arcmin);
        stimset.pedestalSizePix = round(stimset.pedestalSizeAmin/video.pix2arcmin);
        stimset.pedestalLuminance = round(stimset.maxLuminance*stimset.pedestalLum/100);
        stimset.backgroundLuminance = round(stimset.maxLuminance*stimset.bgrLum/100);
        stimset.probeLuminanceMax = round(stimset.pedestalLuminance*(1 + stimset.probeLuminance/100));
        
        frames = zeros( video.minDim,  video.minDim, 1, timing.nUniqueFrames);
        for s = 1:timing.nCoreSteps
            frames(:, :, 1, (s - 1)*timing.nFramesPerStep +1:s*timing.nFramesPerStep) = mkSweepFrames(stimset, timing, video, s);
        end
        
        %% add ciolor dim
        images = uint8(frames);
        
        %% make imageSequence
        imSeq = [];
        updateEveryFrame = timing.nFramesPerStep;
        totalRepeats = timing.framesPerStep/updateEveryFrame;
        for s = 1:timing.nCoreSteps
            activeFrames = (1 + timing.nFramesPerStep*(s - 1):s*timing.nFramesPerStep)';
            imSeq = cat(1, imSeq, repmat(activeFrames, [totalRepeats 1]));
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

    function frames = mkSweepFrames(stimset, timing, video, nStep)
        % get sweep profile
        switch stimset.sweepType
            case 'Pedestal Size'
                % pedestal grows, probe flickers
                pedestalSizePx = getModProfile(stimset.pedestalSizePix(nStep), ...
                    stimset.pedestalSizePix(nStep), timing.nFramesPerStep, stimset.modType);
                probeLum = getModProfile(stimset.probeLuminanceMax(nStep), stimset.pedestalLuminance, ...
                    timing.nFramesPerStep, stimset.modType);
            case 'Probe Contrast'
                probeLum = getModProfile(stimset.probeLuminanceMax(nStep), stimset.pedestalLuminance, ...
                    timing.nFramesPerStep, stimset.modType);
                pedestalSizePx = stimset.pedestalSizePix(nStep)*ones(1, timing.nFramesPerStep);
            case 'Fixed'
                probeLum = getModProfile(stimset.probeLuminanceMax(nStep), stimset.pedestalLuminance, ...
                    timing.nFramesPerStep, stimset.modType);
                pedestalSizePx = stimset.pedestalSizePix(nStep)*ones(1, timing.nFramesPerStep);
                
        end
        
        frames = zeros(video.minDim, video.minDim, timing.nFramesPerStep);
        for nf = 1:timing.nFramesPerStep
            baseElem = makeSquareElement(stimset.baseElemSizePx, stimset.backgroundLuminance, ...
                pedestalSizePx(nf), stimset.pedestalLuminance, ...
                stimset.probeSizePix, probeLum(nf));
            frames(:, :, nf) = repmat(baseElem, stimset.nBaseElems);
        end
        
    end

    function nf = nUniqueFramesPerStep(stimset, timing)
        switch(stimset.modType)
            case 'None'
                nf = timing.updateFramesPerCycle;
            case 'Square'
                nf = timing.updateFramesPerCycle;
            case 'Sawtooth'
                nf = timing.updateFramesPerCycle;
        end
        
    end
    function out = getModProfile(firstVal, lastVal,  nPoints, modProfile)
        switch modProfile
            case 'None'
                out = firstVal*ones(nPoints, 1);
            case 'Square'
                out = [firstVal*ones(round(.5*nPoints), 1); lastVal*ones(nPoints - round(.5*nPoints), 1)];
            case 'Sawtooth'
                out = linspace(firstVal, lastVal, nPoints);
        end
    end
    %% Make Base Element
    function template = makeSquareElement(baseSize, baseLum, pedestalSize, pedestalLum, probeSize, probeLum)
        
        % do a 1/4 template and mirror-flip twice
        probe_first_half = round(probeSize*.5);
        probe_second_half = probeSize - probe_first_half;
        
        pedestal_first_half =  round(pedestalSize*.5);
        pedestal_second_half = pedestalSize - pedestal_first_half;
        
        base_first_half = round(baseSize*.5);
        base_second_half = baseSize - base_first_half;
        
        % first q
        
        template_q1 = baseLum*ones(base_first_half, base_first_half);
        template_q1(1:pedestal_first_half, 1:pedestal_first_half) = pedestalLum;
        template_q1(1:probe_first_half, 1:probe_first_half) = probeLum;
        
        % second q
        template_q2 = baseLum*ones(base_first_half, base_second_half);
        template_q2(1:pedestal_first_half, 1:pedestal_second_half) = pedestalLum;
        if (probe_second_half)
            template_q2(1:probe_first_half, 1:probe_second_half) = probeLum;
        end
        
        % third q
        template_q3 = baseLum*ones(base_second_half, base_first_half);
        template_q3(1:pedestal_second_half, 1:pedestal_first_half) = pedestalLum;
        
        if (probe_second_half)
            template_q3(1:probe_second_half, 1:probe_first_half) = probeLum;
        end
        % q4
        template_q4 = baseLum*ones(base_second_half, base_second_half);
        template_q4(1:pedestal_second_half, 1:pedestal_second_half) = pedestalLum;
        
        if (probe_second_half)
            template_q4(1:probe_second_half, 1:probe_second_half) = probeLum;
        end
        template = [template_q1(end:-1:1, end:-1:1), template_q2(end:-1:1, :); template_q3(:, end:-1:1), template_q4];
    end
end
