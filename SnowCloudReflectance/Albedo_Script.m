% script to compute albedo values from 30 to 2000 um and cosines of solar Z
% from 0 to 1, using SMARTS input for MLW atmosphere at 3 km elevation 

rad = linspace(sqrt(30),sqrt(2000),101);
rad = round(rad.^2);
cZ=0:.02:1;

for c=1:length(cZ)
    atmosP = defaultSMARTSinput('mlw','cosZ',max(cZ(c),.01),'altit',3);
    S = SMARTS295Main(getSMARTShome,'',atmosP);
    T = S.spectralTbl;
    thisTbl = SnowCloudIntgRefl(T.waveL,'nm',[T.HorzDirect T.HorzDiffuse],...
        'snow','cosZ',cZ(c),'radius',rad,'wavelength',[280 4000],...
        'waveu','nm');
    if c==1
        snowTbl = thisTbl;
    else
        snowTbl = [snowTbl; thisTbl];
    end
end    