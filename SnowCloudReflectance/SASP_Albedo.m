function resTbl=SASP_Albedo(whichYear)
% SASP script for albedo values

assert(whichYear>=2005 && whichYear<=2014,...
     'only years 2005 to 2014 are available')
% w = [400 2700; 700 2700]; % wavelengths of clear and red Eppley

w = [305 2800; 780 2800]; %wavelengths of Kipp and Zonen CM21
% variables for this year
L=load('C:\raid\data\nbair\2018albedo_paper\sbbLocation.mat');
data_dir='C:\Users\nbair\Box Sync\SBB\SBB_data\Corrected_Radiation\SASP';
fname=fullfile(data_dir,'SASP_ClearSky_2005-17.csv');
T=readtable(fname);

d = datenum(T.DOY+(datenum(...
    [T.Year ones(height(T),2)])-1))+...
    table2array(T(:,3))/2400;

t=year(d) == whichYear;
fname=fullfile(data_dir,...
    sprintf('sasp_corrected_radiation_%i_fix.dat',...
    whichYear));
T2=readtable(fname);
d2=datenum([table2array(T2(:,1:4)) zeros(height(T2),2)]);

sd_dir=...
    '/Users/nbair/Box Sync/SBB/SBB_data/Precip_SnowDepth/SASP';
fname=fullfile(sd_dir,...
    sprintf('sasp_%i_snowdepth.dat',...
    whichYear));

T3=readtable(fname);
d3=datenum([table2array(T3(:,1:4)) zeros(height(T3),2)]);
sd=table2array(T3(:,6));

%throw out low snow values
t=sd>=0.30;

%clear sky datevals w/ snow depth > 0.30cm
di=intersect(d3(t),d);

%radiation dateval indices for those clear days w/ >= 0.30 cm sd
[d,ix,ixx]=intersect(d2,di);

ac = table2array(T2(ix,9));
ar = table2array(T2(ix,13));

depth = sd(ixx);

[ declin, radiusvector, omega ] = Ephemeris(d+6/24);
mu0 = sunang(L.location.sasp.lat,L.location.sasp.lon, declin, omega);
mu=mu0;

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