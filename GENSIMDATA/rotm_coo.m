function RotM=rotm_coo(Type,EquinoxJD);
%--------------------------------------------------------
% rotm_coo function    Generate rotation matrix for
%                    converting J2000.0 coordinates to
%                    galactic l & b coordinates.
% Input  : - Type of rotation matrix
%            'e'  - equatorial with mean equinox of J2000.0 
%                   to ecliptic with mean ecliptic and equinox J2000.0
%            'E'  - ecliptic with mean ecliptic and equinox J2000.0
%                   to equatorial with mean equinox of J2000.0 
%            'g'  - galactic to equatorial with mean equinox of J2000.0
%            'G'  - equatorial with mean equinox of J2000.0 to galactic.
%            'p'  - precession matrix from mean equinox
%                   of date to mean equinox of J2000.0.
%            'P'  - precession matrix from mean equinox
%                   J2000.0 to mean equinox of date. 
%            'pd' - precession matrix from true equinox
%                   of date to mean equinox of J2000.0.
%            'Pd' - precession matrix from mean equinox
%                   J2000.0 to true equinox of date. 
%          - Equinox of coordinates (in Julian Day),
%            used only in the case of 'p' | 'P' | 'pd' | 'Pd' | 'ed' | 'Ed'
%            In case of 'E' or 'q' if this parameter is
%            not given it is taken as 2451545.0 (=J2000.0)
% Output : - rotation matrix
% Reference : Ex. Supp. to the Astronomical Almanac.
% Tested : Matlab 5.3
%     By : Eran O. Ofek        August 1999  / June 2000  
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------
RADIAN = 180./pi;
J2000  = 2451545.5;

switch Type
 case {'e'}
    Obl  = obliquity(J2000);
    RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)];
    %RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762];
 case {'E'}
    Obl  = obliquity(J2000);
    RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)].';
    %RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762]';   
 case {'g'}
    RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762]';
 case {'G'}
    RotM = [-0.0548755604 +0.4941094279 -0.8676661490; -0.8734370902 -0.4448296300 -0.1980763734; -0.4838350155 +0.7469822445 +0.4559837762];   
 case {'p','P','pd','Pd'}
    T = (EquinoxJD - 2451545.0)./36525.0;
 
    ZetaA  = 0.6406161.*T + 0.0000839.*T.*T + 0.0000050.*T.*T.*T;
    ZA     = 0.6406161.*T + 0.0003041.*T.*T + 0.0000051.*T.*T.*T;
    ThetaA = 0.5567530.*T - 0.0001185.*T.*T - 0.0000116.*T.*T.*T;
    ZetaA  = ZetaA./RADIAN;
    ZA     = ZA./RADIAN;
    ThetaA = ThetaA./RADIAN;

    RotM = zeros(3,3);
    RotM(1,1) = cos(ZetaA).*cos(ThetaA).*cos(ZA) - sin(ZetaA).*sin(ZA);
    RotM(2,1) = cos(ZetaA).*cos(ThetaA).*sin(ZA) + sin(ZetaA).*cos(ZA);
    RotM(3,1) = cos(ZetaA).*sin(ThetaA);
    RotM(1,2) =-sin(ZetaA).*cos(ThetaA).*cos(ZA) - cos(ZetaA).*sin(ZA);
    RotM(2,2) =-sin(ZetaA).*cos(ThetaA).*sin(ZA) + cos(ZetaA).*cos(ZA);
    RotM(3,2) =-sin(ZetaA).*sin(ThetaA);
    RotM(1,3) =-sin(ThetaA).*cos(ZA);
    RotM(2,3) =-sin(ThetaA).*sin(ZA);
    RotM(3,3) = cos(ThetaA);
   
    switch Type
     case {'p'}
        RotM = RotM';
     case {'P'}
        RotM = RotM;
     case {'pd'}
        % calculate nutation matrix
        [Nut, NutMat]=nutation(EquinoxJD);
        RotM = [NutMat']*[RotM'];
     case {'Pd'}
        % calculate nutation matrix
        [Nut, NutMat]=nutation(EquinoxJD);
        RotM = NutMat*RotM;
     otherwise
        error('Unknown rotation matrix type');
    end    
        
        
 otherwise
    error('Unknown rotation matrix type');
end
