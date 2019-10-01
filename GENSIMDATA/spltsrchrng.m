function splitMinMax = spltsrchrng(srchRng,varargin)
%Split a search range matrix along a given parameter
%S = SPLTSRCHRNG(X,P,N)
%Takes a matrix X of maximum and minimum values specifying the search ranges
%for a bunch of parameters (one parameter per row) and splits it into
%matrices where a given parameter's range is split into a specified number of
%uniform bands. The name of the parameter whose range is to be split is
%given by the string P and the number of bands is N. The possible values of
%P are:
%     'SkyPosAlpha','SkyPosDelta', 'AngularFrequency', 'InitialPhase',
%     'LogAmplitude', 'Inclination', 'Polarization'
%
%