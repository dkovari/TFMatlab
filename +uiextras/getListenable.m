function [observable,nonObservable] = getListenable(hndl)
% Return an alphabetical list of Observable properties in handle HNDL
%
% Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 11/11/2015
%
%%% EXAMPLE Property Listener:
% f = figure;
% d = dir;
% a = addlistener(f,?
% ?ApplicationData?,?PostSet?,?
% @(~,~) disp(?appdatachanged?));
%%% This triggers the callback:
% f.ApplicationData.test = d;
%%% This does not. WHY?
% setappdata(f,?test?,d)
 
% Copyright 2015 The MathWorks, Inc.
 
tmp = metaclass(hndl);
properties = tmp.PropertyList;
names = {properties(:).Name}';
observability = {properties(:).SetObservable}';
observability = [observability{:}]';
observable = sort(names(observability));
if nargout > 1
nonObservable = sort(names(~observability));
end