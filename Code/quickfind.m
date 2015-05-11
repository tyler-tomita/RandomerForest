function [loc found]=quickfind(el,list)
% loc=quickfind(el,list)
%
% Finds the index of an element el in a SORTED list really quick.
% loc is the location, found is a boolean indicating whether an exact match
% was found.
%
% It doesn't check if the list is sorted, because that would take O(N) time
% and this is meant to run in O(log(N)).  If you give it an unsorted list,
% it will either return a zero, and error, or the wrong answer.  It is
% guaranteed to stop though, which is nice.
%
% If you give it an element that's not in the sorted list, it interpolates,
% and sets found to false.
% Special cases: if it's smaller than the first element it returns zero.
% If it's bigger then the last element it returns length(list)+1, and of
% course found will be false.
%
% Example usage
% r=sort(rand(1,10000000));
% tic, quickfind(r(7654321),r), toc
%
% Have fun..
%
% Peter
% oconnorp -at- ethz -dot- ch



% If not in bounds, enforce the "negative nearest lower index" rule
if el<list(1), loc=0; found=false; return;
elseif el>list(end), loc=numel(list)+1; found=false; return;
end


n=100; 
% n should be at least 3 to work  Test show it's more or less the same from 
% 30 to 500, then it starts sucking for bigger numbers.
    
st=1;
en=numel(list);
while true
    
    ix=ceil(linspace(st,en,n));
    
    i=find(el<=list(ix),1);
    
    if el==list(ix(i))
        loc=ix(i);
        found=true;
        return;
    else
        st=ix(i-1); en=ix(i);        
        if (en-st)<2 && st~=el,
            % It's not in the list, or the list ain't sorted
            loc=(el-list(st))/(list(en)-list(st))+st;
            found=false;
            return;
        end            
    end    
end
end