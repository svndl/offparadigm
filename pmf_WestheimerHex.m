function pmf_WestheimerHex( varargin )
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
        'View Dist (cm)'	57.3	        'double'  {}
        'Mean Lum (cd)'     11.0	        'double'  {} % under default calibration, this is (0,0,0)
        'Fix Point'         'None'	        'nominal' { 'None' 'Cross' }
        %'Sweep Type'        'Pedestal Size' 'nominal' { 'Fixed' 'Pedestal Size' 'Probe Contrast' }
        'Sweep Type'        'Fixed'         'nominal' { 'Fixed' 'Probe Contrast' }        
        'Step Type'         'Lin Stair'	    'nominal' { 'Lin Stair' 'Log Stair' }
        'Sweep Start'       0.0	            'double'  {}
        'Sweep End'         1.0             'double'  {}
        'Modulation'        'None'          'nominal' { 'None', 'Sawtooth-on', 'Sawtooth-off', 'Square' } %
        }
        %
        
        % "Base" part parameters - paradigm specific parameters that apply to unmodulated parts of the stimulus
        {
        'ModInfo'                   0.0       'integer' {}
        'Bgr lum (% max)'           10        'integer' {} 
        'Pedestal lum (% max)'      50        'integer' {} 
        'Base Element Size (amin)'  50.0      'double'  {}
        'Pedestal Size (% base)'    80        'double'  {}                
		'Magnification Factor'      2.5       'double'	{}      
        }
        
        % "Part1" - parameters that apply to part of stimulus that carries first frequency tag.
        % "Cycle Frames" must be first parameter
        {
        'Cycle Frames'            2.0         'integer'	{}  % framerate(Hz)/stimFreq(Hz) Modulation frequency every other frame (2 Hz)
        'Contrast (pct)'	      20.0        'double'	{}
        'Probe Polarity'         'positive'   'nominal'  {'positive', 'negative'}
        'Probe size (% base)'     20          'double'   {}
        'Sector'                 'All'        'nominal'  {'All', 'f1', 'p1', 'p2'}
        'Visual Field H'         'All'        'nominal'  {'All', 'Upper', 'Lower'}
        'Visual Field V'         'All'        'nominal'  {'All', 'Left', 'Right'}
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
        %'Pedestal Size'    'constant' { { 'B' 'Pedestal Size (amin)' } }
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
         'Sawtooth-on'      'ModInfo'          0.0
         'Sawtooth-off'     'ModInfo'          0.0         
         'Square'           'ModInfo'          0.0
         }
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
        minDim = width_pix;
        
        width_cm = VMVal('imageWidthCm');
        viewDistCm = PVal('S','View Dist (cm)');
        width_deg = 2 * atand( (width_cm/2)/viewDistCm );
        pix2arcmin = ( width_deg * 60 ) / width_pix;
        
        %% Standard output routine
        
        validationMessages = {};
        parametersValid = isempty( validationMessages );
        output = { parametersValid, parameters, validationMessages };
        assignin( 'base', 'output', output );
        
        %% Correct/Message functions
        
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
        stimset.bgrLum = PVal('B', 'Bgr lum (% max)' );
        
        stimset.pedestalLum = PVal('B', 'Pedestal lum (% max)');
        stimset.mfactor = PVal('B', 'Magnification Factor');      
        
        stimset.baseSizeAmin = PVal('B', 'Base Element Size (amin)');
        stimset.pedestalSizePct = GetParamArray('B', 'Pedestal Size (% base)');
        stimset.probeSizePct = PVal('1', 'Probe size (% base)');
        
        stimset.probeLuminance = GetParamArray('1', 'Contrast (pct)');
        stimset.probePolarity =  GetParamArray('1', 'Probe Polarity');
        stimset.probeTable = PVal('1', 'Sector');
        
        try
            stimset.UpDown = PVal('1', 'Visual Field H');    
            stimset.LeftRight = PVal('1', 'Visual Field V');
        catch
            stimset.UpDown ='All';    
            stimset.LeftRight = 'All';            
        end
        save('hexWestheimerInput.mat', 'stimset', 'video', 'stimsetTiming');
        
        [rImSeq, rIms] = westheimerHex(stimset, video, stimsetTiming);
        save('hexWestheimerOutput.mat', 'rImSeq', 'rIms');
                
        isSuccess = true;
        output = { isSuccess, rIms, cast( rImSeq, 'int32') }; % "Images" (single) and "Image Sequence" (Int32)
        %clear rIms
        assignin( 'base', 'output', output )
    end
end

 
