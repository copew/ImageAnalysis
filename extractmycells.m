function C =extractmycells(C)
if iscell(C)
    C = cellfun(@extractmycells, C, 'UniformOutput', 0);
    C = cat(2,C{:});
elseif isempty(C)
    return
else
    C = {C};
end
end