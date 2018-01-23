function [detrenddata]=op_detrend(data,goodi,goodj)
%
% [detrenddata]=op_detrend(data,goodi,goodj);
%
% Detrend optical data. This function is built around the Matlab
% function detrend.m
%
% sister function is op_norm for amplitude normalize 
%
% see also op_norm function
%
% 2003, MWKay
%
detrenddata=data;

if nargin<3
  for i=1:size(data,1)
    for j=1:size(data,2)
      sitedat=reshape(data(i,j,:),1,size(data,3));
      detrenddata(i,j,:)=detrend(sitedat,'linear');
    end
  end
else
  for i=1:length(goodi)
    sitedat=reshape(data(goodi(i),goodj(i),:),1,size(data,3));
    detrenddata(goodi(i),goodj(i),:)=detrend(sitedat,'linear');
  end
end
