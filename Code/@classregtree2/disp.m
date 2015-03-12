function disp(t)
%DISP Display a CLASSREGTREE object.
%   DISP(T) prints the CLASSREGTREE object T.
%
%   See also CLASSREGTREE, CLASSREGTREE/VIEW.

%   Copyright 2006-2011 The MathWorks, Inc.


isLoose = strcmp(get(0,'FormatSpacing'),'loose');
if (isLoose), fprintf('\n'); end

% Get some information about the whole tree
maxnode = numel(t.node);
nd = 1 + floor(log10(maxnode)); % number of digits for node number
varnames = names(t);
if isempty(varnames)
    numCell = textscan(sprintf('%d\n',1:t.npred),'%s\n');
    varnames = strcat('x',numCell{1});
end
isregression = isequal(t.method,'regression');
if isregression
    fprintf(getString(message('stats:classregtree:disp:DecisionTreeForRegression')));
else
    fprintf(getString(message('stats:classregtree:disp:DecisionTreeForClassification')));
end

% Display information about each node
for j=1:maxnode
    if any(t.children(j,:))
        % branch node
        vnum = t.var(j);
        vname = varnames{abs(vnum)};
        cut = t.cut{j};
        kids = t.children(j,:);
        if     strcmp(type(t),'regression')
            Yfit = t.class(j);
            Yfit = num2str(Yfit,'%g');
        elseif strcmp(type(t),'classification')
            Yfit = t.classname{t.class(j)};
        end
        if vnum>0        % continuous predictor "<" condition
            condleft = sprintf('%s<%g',vname,cut);
            condright = sprintf('%s>=%g',vname,cut);
            fprintf('%*d  %s\n',nd,j,getString(message('stats:classregtree:disp:TreeBranch',...
                condleft,kids(1),condright,kids(2),Yfit)));
        else             % categorical predictor, membership condition
            cats = cut{1};
            if isscalar(cats)
                condleft = sprintf('%s=%g',vname,cats);
            else
                set = deblank(num2str(cats,'%g '));
                condleft = sprintf('%s %s {%s}',vname,getString(message('stats:classregtree:disp:ElementInSet')),set);
            end
            cats = cut{2};
            if isscalar(cats)
                condright = sprintf('%s=%g',vname,cats);
            else
                set = deblank(num2str(cats,'%g '));
                condright = sprintf('%s %s {%s}',vname,getString(message('stats:classregtree:disp:ElementInSet')),set);
            end
            fprintf('%*d  %s\n',nd,j,getString(message('stats:classregtree:disp:TreeBranch',...
                condleft,kids(1),condright,kids(2),Yfit)));
        end
    else
        % terminal node, display fit (regression) or class assignment
        if isregression
            fprintf(sprintf('%s  %s %s\n','%*d',getString(message('stats:classregtree:disp:FittedResponse')),'%g'),nd,j,t.class(j));
        else
            fprintf(sprintf('%s  %s %s\n','%*d',getString(message('stats:classregtree:disp:PredictedClass')),'%s'),nd,j,t.classname{t.class(j)});
        end
    end
end
if (isLoose), fprintf('\n'); end
