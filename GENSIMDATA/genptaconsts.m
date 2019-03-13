function constVal = genptaconsts(constName)

%% ======== constants ===========
pc2ly=3.261563777;  % 1 pc=3.26 ly (Julian)
dy2yr=1.0/365.25;  % 1 day=365.25 yr (Julian)  ??
kilo=1.0*10^3;  % kilo 1000

switch constName
    case 'pc21y'
        constVal = 3.261563777;
    case 'dy2yr'
        constVal = 1.0/365.25;
    case 'kilo'
        constVal = 1.0*10^3;
    otherwise
        error('Unrecognized constant name');
end
