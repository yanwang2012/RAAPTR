function Res=cosined(Vec);
%--------------------------------------------------------------------
% cosined function     cosine direction transformation
%                    convert lon & lat to cosine direction
%                    and visa versa. UNITS are in radians!
% input  : - column matrix, with 2 or 3 colums.
%            if 2 colums are given then, the first column
%            is longitude and the second is latitude.
%            and the cosine direction are returned.
%            if three columns are given, then they suppose to
%            be the cosine direction, and the longitude
%            and latitude are returned.
% output : - cosine direction or longitude and latitude.
%    By  Eran O. Ofek           July 1999
%--------------------------------------------------------------------
if (length(Vec(1,:))==2),
   Alpha = Vec(:,1);
   Delta = Vec(:,2);
   Res          = zeros(length(Vec(:,1)),3);
   Res(:,1)     = cos(Alpha).*cos(Delta);
   Res(:,2)     = sin(Alpha).*cos(Delta);
   Res(:,3)     = sin(Delta);
elseif (length(Vec(1,:))==3),
   L1           = Vec(:,1);
   L2           = Vec(:,2);
   L3           = Vec(:,3);
   Res          = zeros(length(Vec(:,1)),2);
   Res(:,1)     = atan2(L2,L1);                 % Alpha
   SLL          = sqrt(L1.^2+L2.^2);
   I0  = find(SLL==0);
   In0 = find(SLL~=0);
   Res(In0,2)     = atan(L3(In0)./SLL(In0));  % Delta
   Res(I0,2)      = sign(L3(I0)).*pi./2;  % Delta
else
   error('only 2/3 columns are allowed');
end
