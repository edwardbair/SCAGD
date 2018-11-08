function [oStruct,varargout] = invertSnowCloudRefl(reflectance,unknowns,varargin )
% [oStruct] = invertSnowCloudRefl(reflectance,unknowns,prescription,prop/val )
% [oStruct] = invertSnowCloudRefl(reflectance,unknowns,substance,prescription,prop/val )
% [oStruct,stats] = invertSnowCloudRefl(reflectance,unknowns,prop/val)
% [ostruct,stats,P] = invertSnowCloudRefl(reflectance,unknowns,prop/val)
%
%solves for snow or cloud properties based on spectral or band-averaged
%   reflectance (i.e. the inverse of SnowCloudSpectralRefl or SnowCloudIntgRefl)
%
%Input
%   reflectance - measured reflectance at wavelengths or band passes
%       corresponding to the input prescription or prop/val pairs
%   unknowns - cell vector of snow or cloud properties to solve for, with
%       possibilities depending on the input substance and the number of
%       reflectance values available (reasonable abbreviations work)
%       for 'snow' - 'radius' or 'ssa', 'fSCA', 'wetness', 'waterEquivalent', 'dust',
%           'dustRadius', 'soot', 'sootRadius', 'corrFactor' (cosS * muS)
%       for 'iceCloud' - 'radius', 'waterEquivalent', 'dust', 'dustRadius', 'soot',
%           'sootRadius'
%       for 'waterCloud' - 'radius', 'waterEquivalent', 'dust', 'dustRadius', 'soot',
%           'sootRadius'
%       for 'mixedCloud' - 'radius', 'waterRadius', 'wetness', 'waterEquivalent',
%           'dust','dustRadius', 'soot', 'sootRadius'
%
%optional arguments all name-value pairs, except substance and prescription
%arguments, if present, must come before the name-value pairs if either or
%both are specified
%   'substance', either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%       (but any unambiguous abbreviation beginning with first letter works)
%   prescription, to use and possibly modify an existing prescription
%Other inputs, name-value pairs, see SetSnowCloud for all possibilities
%for example
%   'cosZ' - cosine of illumination angle on flat surface, scalar, vector, or matrix
%   'muS' - cosine of illumination angle on slope
%   'cosS' - cosine of slope
%   'radius' - effective optical radius of scatterer, scalar, vector, or matrix
%       (if not scalars, 'cosZ' and 'radius' must be same size as wavelength)
%   'ssa' - surface-specific area, kg/m^2, as an alternative to 'radius'
%   'sensor' - any sensor known to SensorTable.m
%   'bands' - spectral bands to consider, can be numeric, character, string,
%       or categorical (default all bands for that sensor)
%   or, instead of 'bands', 'wavelengths' - Nx2 matrix of band passes to
%       integrate across, or vector or matrix of wavelengths for a spectrometer
%   'solutionMethod' - to solve the inversion, either 'lsqnonlin' or
%       'spectralAngle'
%
%Output
%   oStruct - snow or cloud properties, depending on inputs
%Optional output
%   stats - statistics about solution
%   P - prescription at last call to snow reflectance function

narginchk(3,Inf)
nargoutchk(0,3)
passP = struct([]);

% parse inputs (pass to SetSnowCloud)
fscript = SetSnowCloud(varargin{:});
assert(length(reflectance)==length(fscript.wavelength),...
    'length of reflectance vector must equal length of wavelength vector/matrix')
reflectance = double(reflectance);

% parse unknowns
assert(length(unknowns)<=length(reflectance),...
    'length of reflectance vector must be >= number of unknowns')
validStrings = {'fSCA','wetness','dust','dustRadius',...
    'soot','sootRadius','corrFactor','waterRadius','radius',...
    'waterEquivalent','fractionalCoverage','ssa'};
equivStrings = {'fSCA','fractionalCoverage'};

solveFor = cell(length(unknowns),1);
for k=1:length(unknowns)
    matchedstr = validatestring(unknowns{k},validStrings);
    if contains(matchedstr,equivStrings(:),'IgnoreCase',true)
        if contains(matchedstr,equivStrings(:),'IgnoreCase',true)
            solveFor{k} = equivStrings{1};
        end
    else
        solveFor{k} = matchedstr;
    end
end

% make sure wavelenghts long enough to get size of scatterer, if asked for
if any(contains(solveFor,'iceRadius','IgnoreCase',true)) ||...
        any(contains(solveFor,'ssa','IgnoreCase',true)) ||...
        any(contains(solveFor,'waterRadius','IgnoreCase',true))
    assert(convertLengthUnits(max(fscript.wavelength(:)),fscript.waveUnit,'nm')>=1060,...
        'maximum wavelength must be >= %f %s to retrieve size of snow or cloud scatterer',...
        convertLengthUnits(1060,'nm',fscript.waveUnit),fscript.waveUnit)
end
% make sure wavelengths short enough to get dust or soot characterization
if any(contains(solveFor,'dust','IgnoreCase',true)) ||...
        any(contains(solveFor,'soot','IgnoreCase',true))
    assert(convertLengthUnits(min(fscript.wavelength(:)),fscript.waveUnit,'nm')<=700,...
        'minimum wavelength must be <= %f %s to retrieve dust or soot properties',...
        convertLengthUnits(700,'nm',fscript.waveUnit),fscript.waveUnit)
end

% intial values and limits
[x0,lb,ub,extraEndMember] = setBounds(solveFor,fscript);

% solving method depends on input
if fscript.spectrometer
    if size(fscript.R0,2)==2
        R0 = mean(fscript.R0,2);
    else
        R0 = fscript.R0;
    end
    if any(contains(solveFor,'soot'))
        contam = 'soot';
    else
        contam = 'dust';
    end
    passWeight = spectralWeight(fscript.cosZ,fscript.wavelength,...
        fscript.waveUnit,R0,contam);
else
    passWeight = 1;
end

switch fscript.solutionMethod
    % inversion method lsqnonlin uses the signed differences between measurement and model
    case 'lsqnonlin'
        [x,resnorm,residual,exitflag,output,lambda,jacobian] =...
            lsqnonlin(@SnowCloudDiff,x0,lb,ub);
        stats.resnorm = resnorm;
        stats.residual = residual;
        stats.exitflag = exitflag;
        stats.output = output;
        stats.lambda = lambda;
        stats.jacobian = jacobian;
    case 'spectralAngle'
    otherwise
        error('''solutionMethod'' ''%s'' not recognized',fscript.solutionMethod)
end

% put solution into output structure
for k=1:length(solveFor)
    oStruct.(solveFor{k}) = x(k);
    if strcmpi(solveFor{k},'fSCA') && size(fscript.R0,2)>1
        oStruct.otherEndMem = [x(end) 1-x(end)-x(k)];
    end
end

if nargout>1
    varargout{1} = stats;
    if nargout>2
        varargout{2} = passP;
    end
end

    function diffR = SnowCloudDiff(x)
        % difference between snow and cloud reflectance
        
        argc = cell(2*length(solveFor),1);
        v = 1;
        for m=1:length(solveFor)
            if strcmpi(solveFor{m},'fSCA')
                argc{v} = 'fractionalCoverage';
                if extraEndMember
                    argc{v+1} = [x(m) x(end) 1-x(m)-x(end)];
                else
                    argc{v+1} = [x(m) 1-x(m)];
                end
            else
                argc{v} = solveFor{m};
                argc{v+1} = x(m);
            end
            v = v+2;
        end
        % correct for topography either by recomputing effective modeled
        % reflectance (i.e. when viewed from above assuming flat surface)
        % or by recomputing reflectance on slope
        if any(contains(solveFor,'corrFactor','IgnoreCase',true))
            v = contains(solveFor,'corrFactor','IgnoreCase',true);
            m = find(v);
            cf = argc{2*m};
            doCorrFactor = true;
            % correct R0 (effective) for slope and illumination
            if isfield(fscript,'R0')
                if any(fscript.R0~=0)
                    argc{end+1} = 'R0';
                    newR0 = zeros(size(fscript.R0));
                    for c=1:size(fscript.R0,2)
                        newR0(:,c) = fscript.R0(:,c).*fscript.cosZ./cf;
                    end
                    argc{end+1} = newR0;
                end
            end
        else
            doCorrFactor = false;
        end
        if fscript.spectrometer
            [R,passP] = SnowCloudSpectralRefl(fscript,argc{:});
            modelRefl = R.refl;
        else
            [T,passP] = SnowCloudIntgRefl(fscript,argc{:});
            modelRefl = T.reflectance;
        end
        % correct for topography by recomputing effective modeled reflectance (i.e.
        % when viewed from above assuming flat surface)
        if doCorrFactor
            modelRefl = modelRefl.*cf./fscript.cosZ;
        end
        diffR = passWeight.*(reflectance-modelRefl);
    end

end