function [ T] = SnowCloudIntgRefl(irradWavelength,waveUnits,irradiance,varargin)
% [ T] = SnowCloudIntgRefl(irradWavelength,waveUnits,irradiance,substance,[prescription or prop/value pairs\)
%band-integrated reflectance and transmittance of snow or cloud, coupling
%spectral models of snow/cloud reflectance and incoming solar radiation
%in the calculation of snow reflectance across specified wavelenth bands
%
%Required inputs, for the spectrum of the incoming irradiance
%   irradWavelength - vector of wavelengths for the irradiance
%   waveUnits - typically 'mum' or 'nm'
%   irradiance - either vector of same length L as vector of wavelengths,
%       or matrix of size Lx2 where column 1 is the direct irradiance and
%       column 2 is the diffuse irradiance
%
%Variable inputs
%See SetSnowCloud.m for the variables following irradiance to specify the
%snow or cloud properties.
%
%The first values must be the substance, either 'snow', 'iceCloud',
%'waterCloud', or 'mixedCloud' (any unambiguous abbreviation beginning
%with first letter works) or you can specify the prescription output from
%SetSnowCloud.
% Wavelength properties, either 'wavelength' or 'sensor'/'band' must be
%   specified, but not both
% 'waveUnit' – units for wavelength, default 'mum', applies whether
%   'wavelength' or 'sensor'/'band' are set
% 'wavelength' – Nx2 matrix for multispectral sensors, or use [.28 4] or
%   [280 4000] to get full spectrum albedo
% 'sensor' – instead of 'wavelength', can specify any multispectral sensor
%   in the SensorTable.m function
% 'bands' – bands of the sensor, either numeric vector, cell vector, or
%   categorical vector of bands,or, if omitted, all bands for that sensor
% 'ignoreSolar' – If false (default), solar radiation accounted for unless
%   outside range of input solar radiation values from SMARTS.
%   If true, ignores solar radiation and just provides band-average
%   reflectivity (this is needed to calculate emissivity around 4 um).
%
%Values for 'radius' and 'cosZ' can be specified as scalars or vectors, but if
%vectors they must be the same size. to get all combinations of radius/cosZ
%values, use meshgrid or ndgrid first and convert the results to vectors.
%Make sure you use the same solar geometry in the software that generates
%the irradiance values.
%
%Output
% T - table of snow or cloud reflectance and transmittance, dimensionless,
%   same height as number of bands x number of radii (which equals number
%   of cosZ)

%%
narginchk(4,Inf)
nargoutchk(0,1)

SnowCloudP = SetSnowCloud(varargin{:});
assert(~SnowCloudP.spectrometer,'code %s not suitable for a spectrometer',mfilename);

%wavelength units must be the same
assert(strcmpi(SnowCloudP.waveUnit,waveUnits) ||...
    (contains(SnowCloudP.waveUnit,'um','IgnoreCase',true) &&...
    contains(waveUnits,'um','IgnoreCase',true)),...
    'wavelength units for the irradiance wavelengths and the snow/cloud wavelengths must be the same')

% generate reflectance spectrum for the snow or cloud
% hold the wavelength ranges for later integration
bandPass = SnowCloudP.wavelength;
minW = min(bandPass(:));
maxW = max(bandPass(:));
if contains(char(SnowCloudP.substance),'water','IgnoreCase',true)
    [~,wv] = RefractiveIndex([],'water',SnowCloudP.waveUnit);
else
    [~,wv] = RefractiveIndex([],'ice',SnowCloudP.waveUnit);
end
k1 = find(wv<=minW,1,'last');
k2 = find(wv>=maxW,1,'first');
wv = wv(k1:k2);
P1 = SnowCloudP;

% n-d grid of radius cosZ and wavelength
snow = categorical({'snow'});
waterCloud = categorical({'waterCloud'});
iceCloud = categorical({'iceCloud'});
mixedCloud = categorical({'mixedCloud'});
switch P1.substance
    case {snow,iceCloud,mixedCloud}
        [r,cz] = ndgrid(unique(P1.iceRadius),unique(P1.cosZ));
        r = r(:);
        cz = cz(:);
        [radius,wave] = ndgrid(r,wv);
        P1.iceRadius = radius;
    case waterCloud
        [r,cz] = ndgrid(unique(P1.iceRadius),unique(P1.cosZ));
        r = r(:);
        cz = cz(:);
        [radius,wave] = ndgrid(r,wv);
        P1.waterRadius = radius;
    otherwise
        error('substance ''%s'' not recognized',char(P1.substance))
end
[P1.cosZ,~] = ndgrid(cz,wv);
P1.wavelength = wave;

% direct reflectance
M = SnowCloudSpectralRefl(P1);
refl = M.refl;
% diffuse reflectance
if size(irradiance,2)==2
    P2 = P1;
    P2.cosZ = [];
    M = SnowCloudSpectralRefl(P2);
    refl = [refl M.refl];
end

reflDir = reshape(refl(:,1),size(P1.wavelength));
if size(refl,2)==2
    reflDif = reshape(refl(:,2),size(P1.wavelength));
end
for r=1:size(P1.cosZ,1)
    cosZ = unique(P1.cosZ(r,:));
    w = unique(P1.wavelength(r,:))';
    if exist('reflDif','var')
        thisRefl = [reflDir(r,:)' reflDif(r,:)'];
    else
        thisRefl = reflDir(r,:)';
    end
    switch P1.substance
        case {snow,iceCloud,mixedCloud}
            rad = unique(P1.iceRadius(r,:));
            radName = 'iceRadius';
        case waterCloud
            rad = unique(P1.waterRadius(r,:));
            radName = 'waterRadius';
    end
    for b=1:size(bandPass,1)
        R = bandPassReflectance(w,P1.waveUnit,thisRefl,irradiance,...
            'bandPass',bandPass(b,:),'irradWavelength',irradWavelength);
        if ~isempty(P1.sensor)
            thisTbl = table({P1.sensor},P1.bands(b),bandPass(b,:),rad,cosZ,R,...
                'VariableNames',...
                {'sensor','band','bandPass',radName,'cosZ','reflectance'});
            thisTbl.Properties.VariableUnits =...
                {'','',waveUnits,P1.sizeUnit,'',''};
        else
            thisTbl = table(bandPass(b,:),rad,cosZ,R,'VariableNames',...
                {'bandPass',radName,'cosZ','reflectance'});
            thisTbl.Properties.VariableUnits =...
                {waveUnits,P1.sizeUnit,'',''};
        end
        if r==1 && b==1
            T = thisTbl;
        else
            T = [T; thisTbl]; %#ok<AGROW>
        end
    end
end
end