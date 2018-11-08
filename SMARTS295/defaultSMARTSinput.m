function prescription = defaultSMARTSinput(RefAtmos,varargin)
%defaultSMARTSinput - generate input for SMARTS295Main with default values
%
%Required input
%   RefAtmos - reference atmosphere to use, options are (case-insensitive):
%       'mlw' - mid-latitude winter
%       'sas' - Sub-Arctic summer
%       others to be added
%Optional input, name/value pairs
%   'cosZ' - cosine of solar zenith
%   'altit' - elevation of surface, km (default 3 for 'mlw', 0 for 'sas')
%   'height' - elevation above surface, km (default 0)
%
%For other options, use the more general routine SetSMARTS295

p = inputParser;
defaultZen = 48.19;
defaultCosZ = cosd(defaultZen);
defaultALTIT = 0;
defaultHEIGHT = 0;
addRequired(p,'RefAtmos',@ischar)
addParameter(p,'cosz',defaultCosZ,...
    @(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1)
addParameter(p,validatestring('altit',{'alt','elev','altit'}),...
    defaultALTIT,@(x) isscalar(x) && isnumeric(x) && x<=70)
addParameter(p,validatestring('height',{'ht','height'}),...
    defaultHEIGHT,@(x) isscalar(x) && isnumeric(x) && x<=70)
parse(p,RefAtmos,varargin{:});

switch lower(p.Results.RefAtmos)
    case 'mlw'
        argc = {'COMNT','mid-latitude winter','ISPR',1,'SPR',700,...
            'ALTIT',3,'LATIT',37.6,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',30,...
            'IPRT',2,'IMASS',0,'ZENIT',defaultZen,'ATMOS','MLW','IOUT',...
            [1 2 3 5 21 30]};
    case 'sas'
        argc = {'COMNT','sub-Arctic summer','ISPR',1,...
            'LATIT',70,'IATMOS',1,'IH2O',1,'ISPCTR',0,'IALBDX',0,...
            'IPRT',2,'IMASS',0,'ZENIT',defaultZen,'ATMOS','SAS','IOUT',...
            [1 2 3 5 21 30]};
    otherwise
        error('RefAtmos ''%s'' unknown',p.Results.RefAtmos)
end
% adjust defaults if specified
cosZ = p.Results.cosz;
if cosZ~=defaultCosZ
    t = strcmpi(argc,'ZENIT');
    if nnz(t)
        k = find(t);
        argc{k+1} = acosd(cosZ);
    else
        argc = cat(2,argc,{'ZENIT',acosd(cosZ)});
    end
end
altit = p.Results.altit;
if altit~=defaultALTIT
    t = strcmpi(argc,'ALTIT');
    if nnz(t)
        k = find(t);
        argc{k+1} = p.Results.altit;
    else
        argc = cat(2,argc,{'ALTIT',p.Results.altit});
    end
end
height = p.Results.height;
if height~=defaultHEIGHT
    t = strcmpi(argc,'HEIGHT');
    if nnz(t)
        k = find(t);
        argc{k+1} = p.Results.height;
    else
        argc = cat(2,argc,{'HEIGHT',p.Results.height});
    end
end

%set the prescription
prescription = SetSMARTS295(argc{:});
end