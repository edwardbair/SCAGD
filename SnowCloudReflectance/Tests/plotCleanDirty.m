function [X] = plotCleanDirty(contam,sensor,sizeUnit,waveUnit)
% [X] = plotCleanDirty(contam,sensor,sizeUnit,waveUnit)
%make a plot of clean vs dirty snow either for spectrometer or
%multispectral sensor
%return values in a structure, depending on inputs

% set up the inputs
S = SnowCloudLimits;
r = linspace(sqrt(convertUnits(50,'um',sizeUnit)),sqrt(convertUnits(1500,'um',sizeUnit)),11).^2;
switch contam
    case 'dust'
        contamR = rand(1)*(S.dustRadius(2)-S.dustRadius(1))+S.dustRadius(1);
        contamR = convertUnits(contamR,S.unitsSize,sizeUnit);
        conc = logspace(log10(2*S.dust(1)),log10(S.dust(2)),5);
    case 'soot'
        contamR = rand(1)*(S.sootRadius(2)-S.sootRadius(1))+S.sootRadius(1);
        contamR = convertUnits(contamR,S.unitsSize,sizeUnit);
        conc = linspace(S.soot(1),S.soot(2),5);
    otherwise
        error('contam ''%s'' not recognized',contam)
end

% clean snow
for k=1:length(r)
    if k==1
        [R,Pc] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'sizeUnit',sizeUnit,'waveUnit',waveUnit,'radius',r(k));
        refl = zeros(length(Pc.wavelength),length(r));
    else
        [R,Pc] = SnowCloudSpectralRefl(Pc,'radius',r(k));
    end
    refl(:,k) = R.refl;
end
cleanR = refl;
% dirty snow
% first with minimum radius, then with max
for k=1:length(conc)
    rname = [contam 'Radius'];
    if k==1
        [R,Pd] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'sizeUnit',sizeUnit,'waveUnit',waveUnit,'radius',min(r),...
            contam,conc(k),rname,contamR);
        refl = zeros(length(Pd.wavelength),length(conc)*2);
    else
        [R,Pd] = SnowCloudSpectralRefl(Pd,contam,conc(k));
    end
    refl(:,k) = R.refl;
end
for k=1:length(conc)
    rname = [contam 'Radius'];
    if k==1
        [R,Pd] = SnowCloudSpectralRefl('snow','cosZ',cosd(50),'sensor',sensor,...
            'sizeUnit',sizeUnit,'waveUnit',waveUnit,'radius',max(r),...
            contam,conc(k),rname,contamR);
    else
        [R,Pd] = SnowCloudSpectralRefl(Pd,contam,conc(k));
    end
    refl(:,k+length(conc)) = R.refl;
end
dirtyR = refl;

% absorption coefficient of ice
[N,w] = RefractiveIndex([],'ice','nm');
k1 = find(w<min(Pd.wavelength),1,'last');
k2 = find(w>max(Pd.wavelength),1,'first');
w = w(k1:k2);
N = N(k1:k2);

% plot
hdirty=plot(Pd.wavelength,dirtyR,'r','LineWidth',1.5);
hold on;
hclean=plot(Pc.wavelength,cleanR,'b','LineWidth',1);
ylabel('spectral albedo')
xlabel(['wavelength, ' waveUnit])
xlim([convertLengthUnits(300,'nm',waveUnit) convertLengthUnits(2600,'nm',waveUnit)])
yyaxis right
habs=plot(w,imag(N),'k','LineWidth',2);
hglines = [hclean(1) hdirty(1) habs];
ylabel('absorption coefficient')
u = '{\mu}m';
% cleanLabel = ['clean snow, r = ' num2str(round(r(1))) ' to ' num2str(round(r(end))) ' ' u];
cleanLabel = sprintf('clean snow, r = %d to %d %s',...
    round(r(1)), round(r(end)),u);
dustyLabel = sprintf('dusty snow, %d to %d ppm at r = %d & %d %s',...
    round(conc(1)*1e6),round(conc(end)*1e6),round(min(r)),round(max(r)),u);
absLabel = 'ice absorption coefficient';
legend(hglines,cleanLabel,dustyLabel,absLabel,'Location','Best')


%output
X.clearRefl = cleanR;
X.dirtyRefl = dirtyR;
X.Pc = Pc;
X.Pd = Pd;
X.radius = r;
X.concentration = conc;
end


