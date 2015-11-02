function [varargout] = subsasgn(varargin)
%SUBSASGN Subscripted reference for a CLASSREGTREE object.
%   Subscript assignment is not allowed for a CLASSREGTREE object.

%   Copyright 2006-2007 The MathWorks, Inc.


error(message('stats:classregtree:subsasgn:NotAllowed'))
