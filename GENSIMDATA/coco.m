function [OutList,TotRot]=coco(InList,InCooType,OutCooType,InUnits,OutUnits)
%--------------------------------------------------------
% coco function    General coordinate convertor.
%                  Convert/precess coordinate from/to
%                  Equatorial/galactic/Ecliptic.
% Input  : - Matrix of input coordinates.
%            First column for Long/RA and second column
%            for lat/Dec.
%          - Type of input coordinates.
%            'j####.#' - equatorial, mean of date (default 'j2000.0'). [Julian].
%            'J####.#' - equatorial, true of date (default 'j2000.0'). [Julian].
%            'g' - J2000.0 galactic.
%            'e' - Ecliptic with J2000.0 equinox.
%          - Type of outpt coordinates.
%            'j####.#' - equatorial, mean of date (default 'j2000.0'). [Julian].
%            'J####.#' - equatorial, true of date (default 'j2000.0'). [Julian].
%            'g' - J2000.0 galactic. (default)
%            'e' - Ecliptic with J2000.0 equinox.
%          - Units for input coordinates.
%            'r' - radians. (default)
%            'd' - degrees.
%            'h' - hours/deg.
%          - Units for outpu coordinates.
%            'r' - radians. (default)
%            'd' - degrees.
%            'h' - hours/deg.
% Output : - Matrix of output coordinates.
%          - Total rotation matrix.
% Tested : Matlab 5.3
%     By : Eran O. Ofek             Febuary 2000 / June 2000 / December 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%--------------------------------------------------------
RADIAN = 180./pi;

if (nargin==1),
   InCooType  = 'j2000.0';
   OutCooType = 'g';
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==2),
   OutCooType = 'g';
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==3),
   InUnits    = 'r';
   OutUnits   = 'r';
elseif (nargin==4),
   OutUnits   = 'r';
elseif (nargin==5),
   % do nothing
else
   error('Illigal number of input arguments');
end

LenInType  = length(InCooType);
LenOutType = length(OutCooType);

if (LenInType>1),
   InEquinox  = str2num(InCooType(2:LenInType));
   InEquinoxJD = 2451545.5 + 365.25.*(InEquinox - 2000);
   InCooType  = InCooType(1);
end

if (LenOutType>1),
   OutEquinox = str2num(OutCooType(2:LenOutType));
   OutEquinoxJD = 2451545.5 + 365.25.*(OutEquinox - 2000);
   OutCooType  = OutCooType(1);
end

InCoo   = zeros(size(InList));
OutCoo  = zeros(size(InList));
OutList = zeros(size(InList));

switch InUnits
 case {'r'}
    InCoo = InList;
 case {'d'}
    % convert deg. to radians
    InCoo = InList./RADIAN;
 case {'h'}
    % convert h/d to radians
    InCoo(:,1) = InList(:,1).*15./RADIAN;
    InCoo(:,2) = InList(:,2)./RADIAN;
 otherwise
    error('Unknown type of input units');
end

% convert coordinates to direction cosines
InCosDir = cosined(InCoo);


RotM1 = diag([1 1 1]);

% calculate the first rotation matrix
switch InCooType
 case {'j','J'}
    if (InEquinox~=2000.0),
       % precess coordinates to J2000.0
       switch InCooType
       case {'j'}
           % mean equinox ...
           RotM1 = rotm_coo('p',InEquinoxJD);
        case {'J'}
           % true equinox ...
           RotM1 = rotm_coo('pd',InEquinoxJD);
        otherwise
           error('Illegal InCooType');
       end
    end
 case {'g'}
    % convert to Equatorial J2000.0
    RotM1 = rotm_coo('G',2451545.5);
 case {'e'}
    % convert to Equatorial J2000.0
    RotM1 = rotm_coo('E',2451545.5);
 otherwise
    error('Unknown input coordinaytes type');
end

RotM2 = diag([1 1 1]);
% calculate the second rotation matrix
switch OutCooType
 case {'j','J'}
    if (OutEquinox~=2000.0),
       % precess coordinates from J2000.0
       switch OutCooType
        case {'j'}   
           % mean equinox ...
           RotM2 = rotm_coo('P',OutEquinoxJD);
        case {'J'}   
           % true equinox ...
           RotM2 = rotm_coo('Pd',OutEquinoxJD);
        otherwise
           error('Illegal OutCooType');
       end           
    end
 case {'g'}
    % convert to galactic
    RotM2 = rotm_coo('g',2451545.5);
 case {'e'}
    % convert to ecliptic
    RotM2 = rotm_coo('e',2451545.5);
 otherwise
    error('Unknown output coordinaytes type');
end

% rotate coordinates
TotRot = RotM2*RotM1;
OutCosDir = TotRot*[InCosDir.'];


% convert coordinates from direction cosines
OutCoo = cosined([OutCosDir.']);


switch OutUnits
 case {'r'}
    OutList = OutCoo;
 case {'d'}
    % convert radians to deg.
    OutList = OutCoo.*RADIAN;
 case {'h'}
    % convert radians to h/d
    OutList(:,1) = OutCoo(:,1).*RADIAN./15;
    OutList(:,2) = OutCoo(:,2).*RADIAN;
 otherwise
    error('Unknown type of output units');
end






