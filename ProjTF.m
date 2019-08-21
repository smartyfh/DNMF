function Dt = ProjTF(D)
% Projection onto the set of Tight Frames
% D = Input Frame
% Dt = Output tight frame

% Constrained Analysis Operator Learning for Cosparse Signal Modelling
% Written by Mehrdad Yaghoobi, version 1.0                                    
%     
% Copyright 2013 Mehrdad Yaghoobi, Sangnam Nam, Remi Gribonval and Mike E. Davies
% 
% For all details please refer to README.m
%
% This software is a free software distributed under the terms of the GNU 
% Public License version 3 (http://www.gnu.org/licenses/gpl.txt). You can 
% redistribute it and/or modify it under the terms of this licence, for 
% personal and non-commercial use and research purpose. 

[U,S,V] = svd(D);
[M,N] = size(D);
Dt = U*eye(size(S))*V';
