function resTbl=CUES_Albedo(A,whichYear)
% resTbl=CUES_Albedo(whichYear)
% CUES script for albedo values
% A = output from 'Tests\2011to2017albedo_raw.mat'

yrNo = whichYear-2010; % i.e. year 1 is 2011
assert(whichYear>=2011 && whichYear<=2017,...
    'only years 2011 to 2017 are available')
w = [400 2700; 700 2700]; % wavelengths of clear and red Eppley

% variables for this year
ac = A.albedo_c{yrNo};
ar = A.albedo_red_c{yrNo};
d = A.matdates{yrNo};
mu = A.mu{yrNo};
mu0 = A.mu0{yrNo};
depth = A.depth{yrNo};

% keep only the values where neither albedo is NaN;
t = isnan(ac) | isnan(ar);
ac = ac(~t);
ar = ar(~t);
d = d(~t);
mu = mu(~t);
mu0 = mu0(~t);
depth = depth(~t);

% estimate grain size, dust, and deltavis for each day
for k=1:length(ac)
    P = defaultSMARTSinput('mlw','cosZ',mu0(k));
    X = SMARTS295Main(getSMARTShome,P);
    Tatm = X.spectralTbl;
    outStruct = invertCUES(mu(k),Tatm.waveL,Tatm.HorzDirect,Tatm.HorzDiffuse,...
        w,[ac(k) ar(k)]);
    %determine clean snow int. reflectance over MODDRFS wl (350 to
    %876 nm) and subtract dirty snow int. refl. to get deltavis
    idx1=find(Tatm.waveL==350);
    idx2=find(Tatm.waveL==876);
    Rclean = SnowCloudIntgRefl(Tatm.waveL(idx1:idx2),'nm',...
        [Tatm.HorzDirect(idx1:idx2),Tatm.HorzDiffuse(idx1:idx2)],...
        'snow','radius',outStruct.RadiusDust.radius,'dust',0,...
        'cosZ',mu0(k),'wavelength',Tatm.waveL([idx1 idx2])','waveU','nm');
    Rdirty = SnowCloudIntgRefl(Tatm.waveL(idx1:idx2),'nm',...
        [Tatm.HorzDirect(idx1:idx2),Tatm.HorzDiffuse(idx1:idx2)],...
        'snow','radius',outStruct.RadiusDust.radius,'dust',...
        outStruct.RadiusDust.dust,...
        'cosZ',mu0(k),'wavelength',Tatm.waveL([idx1 idx2])','waveU','nm');
    deltavis = Rclean.reflectance-Rdirty.reflectance;
    dt = datetime(d(k),'ConvertFrom','datenum');
    dt.Format = 'yyyy-MMM-dd'' ''HH:mm';
    O = outStruct.RadiusDust;
    thisTbl = table(dt,mu0(k),mu(k),...
        [ac(k) ar(k)],O.resnorm,[ac(k)+O.residual(1) ar(k)+O.residual(2)],...
        O.radius,O.dust,deltavis,O.exitflag,depth(k),...
        'VariableNames',{'date','mu0','mu','measAlbedo','resnorm','modelAlbedo',...
        'radius','dust','deltavis','exitflag','depth'});
    if k==1
        resTbl = thisTbl;
    else
        resTbl = [resTbl; thisTbl]; %#ok<AGROW>
    end
end
end