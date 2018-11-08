function [wt,varargout]=spectralWeight(cosZ,lambda,lambdaUnits,backR,contam)
% [wt [,snowDiff,backDiff,solarTrans]]=spectralWeight(cosZ,lambda,lambdaUnits,backR,contam)
% wt=spectralWeight(cosZ,lambda,lambdaUnits,backR,contam)
% weights based on differences across grain size and contaminant amounts for snow
% and between snow reflectance and the background
%
%Input
% cosZ - cosine of illumination angle
% lambda - wavelength vector
% lambdaUnits - typically 'um', 'mum', or 'nm'
% backR - spectral reflectance of snow endmember
% contam - 'dust' or 'soot'
%
%Output
% wt - spectral weights for use in least squares analysis
%Optional output - components of the weight
% snowDiff - normalized difference between fine, clean and coarse, dirty snow
% backDiff - normalized difference between snow and background
% solarTrans - solar radiation transmittance


S = SnowCloudLimits();
narginchk(5,5)
nargoutchk(0,4)

%fine, clean snow
R1 = SnowCloudSpectralRefl('snow','cosZ',cosZ,'radius',S.snowRadius(1),...
    'sizeUnit',S.unitsSize,'wavelength',lambda,'waveUnit',lambdaUnits);
% coarse clean snow
R2 = SnowCloudSpectralRefl('snow','cosZ',cosZ,'radius',S.snowRadius(2),...
    'sizeUnit',S.unitsSize,'wavelength',lambda,'waveUnit',lambdaUnits);
%coarse, dirty or sooty snow
switch contam
    case 'dust'
        conc = S.dust(2)/2;
    case 'soot'
        conc = S.soot(2)/2;
    otherwise
        error('contam argument must be ''dust'' or ''soot''')
end
% coarse dirty snow
R3 = SnowCloudSpectralRefl('snow','cosZ',cosZ,'radius',S.snowRadius(2)/2,...
    'sizeUnit',S.unitsSize,'wavelength',lambda,'waveUnit',lambdaUnits,...
    contam,conc);

% snow vs background
backDiff = normalize(abs(R1.refl-backR)+abs(R2.refl-backR)+abs(R3.refl-backR),'range');
% weights for the spectral absorptions are doubled
snowDiff = normalize((R1.refl-R2.refl)+(R2.refl-R3.refl),'range');

% solar transmittance
solarTrans = SolarScale(lambda,'units',lambdaUnits,'location','ratio');
w = double(normalize(max(snowDiff,backDiff).*solarTrans,'range'));

if isrow(lambda)
    F = fit(lambda',w','pchipinterp');
else
    F = fit(lambda,w,'pchipinterp');
end
wt = F(lambda);
wt = normalize(wt,'range');
if isrow(wt)
    wt = wt';
end

%optional output
for k=2:nargout
    n = k-1;
    switch n
        case 1
            varargout{n} = snowDiff; %#ok<*AGROW>
        case 2
            varargout{n} = backDiff;
        case 3
            varargout{n} = solarTrans;
    end
end
end