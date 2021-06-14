function [permutation,partsize,partbegin,partend] = OneLevelKwayNestedDissection(A,level)

n = size(A, 2);
permutation = zeros(1, n);
nbPartition = 2^level;

partbegin = zeros(1, nbPartition + 1);
partend   = zeros(1, nbPartition + 1);
partsize  = zeros(1, nbPartition + 1);

[perm,~,sizes] = metismex('NodeNDP', A, nbPartition);

% Get the sizes blocks
x = 1;
for j = 1 : level + 1
    for i = 1 : 2^(level + 1 - j)
        psize{j}{i} = sizes(x);
        x = x + 1;
    end
end

% Set the begin indices of blocks
pbegin = idx_process(1, level + 1, 0, psize); % begin_idx of original NDP system

% Put together all the separators
idx = 1;

for i = 1 : nbPartition
   partbegin(i) = idx;
   permutation(idx : idx + psize{1}{i} - 1)=(pbegin{1}{i} : pbegin{1}{i} - 1 + psize{1}{i});
   idx = idx + psize{1}{i};
end

partbegin(nbPartition + 1) = idx;
for j = 2 : level + 1
    for i = 1 : 2^(level + 1 - j)
        permutation(idx : idx + psize{j}{i} - 1) = (pbegin{j}{i} : pbegin{j}{i} - 1 + psize{j}{i});
        idx = idx + psize{j}{i};
    end
end

for i = 1 : nbPartition
   partsize(i) = psize{1}{i};
   partend(i)  = partbegin(i) + partsize(i) - 1;
end

partsize(nbPartition + 1) = 0;
for j = 2 : level + 1
    for i = 1 : 2^(level + 1 - j)
        partsize(nbPartition + 1) = partsize(nbPartition + 1) + psize{j}{i};
    end
    partend(nbPartition + 1) = partbegin(nbPartition + 1) + partsize(nbPartition + 1) - 1;

end

permutation = perm(permutation);
end

function gbegin = idx_process(node, level, set, psizes)

global gbegin;
global gsizes;
global idx;

if set == 0
    idx = 1;
    gsizes = psizes;
end


if level ~= 1
    idx_process(2 * node - 1, level - 1, 1);
    idx_process(2 * node, level - 1, 1);
end

gbegin{level}{node} = idx;
idx = idx + gsizes{level}{node};
end
