clc();
clear();
close('all');
% lon lat vel jid
data = csvread('../data.csv');
data = filt_by_velocity(data);
tree = build(data);
graph = cell(size(data,1),1);
eps = 8;
minpts = 50;
for i = 1:size(data,1)
    result = sort(range(tree,data,data(i,1:2),eps));
    dist = distance(data(i,1),data(i,2),data(result,1),data(result,2));
    j = 1;
    while j <= numel(result)
        data_index = result(j);
        jid = data(data_index,4);
        min_distance = dist(j);
        min_index = result(j);
        k = j+1;
        while k <= numel(result) && data(result(k),4) == jid
            if dist(k) < min_distance || (dist(k) == min_distance && result(k) < min_index)
                min_distance = dist(k);
                min_index = result(k);
            end
            k = k+1;
        end
        if jid == data(i,4)
            min_index = i;
        end
        graph{i}(1,end+1) = min_index;
        j = k;
    end
    if i/5000 == floor(i/5000)
        assert(isequal(result,find(distance(data(i,1),data(i,2),data(:,1),data(:,2)) <= eps)));
        disp([num2str(i),'/',num2str(size(data,1))]);
    end
end
pts = cellfun(@numel,graph);
cores = find(pts>=minpts);
i = qsort(pts(1:numel(cores)));
cores = cores(i);
group = zeros(size(pts));
group_id = 0;
%{
rust = csvread('../rust_tree.csv');
assert(all(tree.nodes(:,1) == rust(:,1)+1));
assert(all(tree.nodes(:,2) == rust(:,2)));
assert(all(tree.nodes(:,3) == rust(:,3)));
assert(all(tree.nodes(:,4) == rust(:,4)));
assert(all(tree.index == rust(:,5)+1));
%}
for idx = 1:numel(cores)
    i = cores(idx);
    if numel(graph{i}) >= minpts && group(i) == 0
        group_id = group_id+1;
        group = dfs(graph, minpts, group, i, group_id);
    end
end
fid = fopen('../ans_matlab.csv','w+');
for i = 1:size(data,1)
    fprintf(fid, '%f,%f,%d\n', data(i,1), data(i,2), group(i));
end
fclose(fid);
function group = dfs(graph, minpts, group, i, group_id)
    if group(i) ~= 0
        return
    end
    group(i) = group_id;
    if numel(graph{i}) >= minpts
        for j = graph{i}
            group = dfs(graph, minpts, group, j, group_id);
        end
    end
end
function dist = distance(lon1,lat1,lon2,lat2)
    sin_dlon_2 = sin((pi ./ 180) .* (lon1 - lon2) ./ 2);
    sin_dlat_2 = sin((pi ./ 180) .* (lat1 - lat2) ./ 2);
    cos_lat1 = cos((pi ./ 180) .* lat1);
    cos_lat2 = cos((pi ./ 180) .* lat2);
    muladd = (sin_dlat_2 .* sin_dlat_2) + (cos_lat1 .* cos_lat2) .* (sin_dlon_2 .* sin_dlon_2);
    dist = 2.*asin(sqrt(muladd)).*6371000;
end
function result = range(tree,data,coord,maxDistance)
    nodes = tree.nodes;
    index = tree.index;
    stack = [1;zeros(17,1)];
    stacklen = 1;
    result = zeros(10,1);
    resultlen = 0;
    while stacklen > 0
        node_index = stack(stacklen);
        stacklen = stacklen-1;
        data_index = index(node_index);
        if distance(data(data_index,1),data(data_index,2),coord(1),coord(2)) <= maxDistance
            if resultlen == numel(result)
                tmp = zeros(2*resultlen,1);
                tmp(1:resultlen) = result;
                result = tmp;
            end
            resultlen = resultlen+1;
            result(resultlen) = data_index;
        end
        k = nodes(node_index,1);
        x = nodes(node_index,2);
        lo = nodes(node_index,3);
        hi = nodes(node_index,4);
        alter = coord(:);
        alter(k) = x;
        if hi ~= 0 && ~(coord(k) <= x && distance(coord(1),coord(2),alter(1),alter(2)) >= maxDistance)
            stacklen = stacklen+1;
            stack(stacklen,1) = hi;
        end
        if lo ~= 0 && ~(coord(k) >= x && distance(coord(1),coord(2),alter(1),alter(2)) >= maxDistance)
            stacklen = stacklen+1;
            stack(stacklen,1) = lo;
        end
    end
    result = result(1:resultlen);
end
function array = nth_element(array, begin, nth, end_, less)
    while end_-begin >= 1
        n = end_-begin;
        p = begin;
        a = begin+floor(n/2);
        b = end_-1;
        if less(array(p), array(a))
            t = array(p);
            array(p) = array(a);
            array(a) = t;
        end
        if less(array(b), array(p))
            t = array(p);
            array(p) = array(b);
            array(b) = t;
        end
        if less(array(p), array(a))
            t = array(p);
            array(p) = array(a);
            array(a) = t;
        end
        a = begin;
        b = end_;
        while a+1 < b
            if less(array(a+1), array(a))
                t = array(a+1);
                array(a+1) = array(a);
                array(a) = t;
                a = a+1;
            else
                t = array(a+1);
                array(a+1) = array(b-1);
                array(b-1) = t;
                b = b-1;
            end
        end
        if nth == a
            break
        end
        if nth < a
            end_ = a;
        else
            begin = b;
        end
    end
end
function var = variance(data, index, begin, end_, k)
    n = end_-begin;
    mean = 0;
    for i = begin:end_-1
        x = data(index(i),k);
        mean = mean+x/n;
    end
    var = 0;
    for i = begin:end_-1
        x = data(index(i),k);
        var = var+(x-mean)*(x-mean)/(n-1);
    end
end
function tree = build(data)
    nodes = NaN(size(data,1),4);
    index = (1:size(data,1)).';
    stack = [1; size(data,1)+1; 0];
    stacklen = 1;
    nodeslen = 0;
    while stacklen > 0
        begin = stack(1,stacklen);
        end_ = stack(2,stacklen);
        parent = stack(3,stacklen);
        stacklen = stacklen-1;
        n = end_-begin;
        if n == 0
            if parent ~= 0
                if isnan(nodes(parent,3))
                    nodes(parent,3) = 0;
                else
                    nodes(parent,4) = 0;
                end
            end
            continue
        end
        var1 = variance(data, index, begin, end_, 1);
        var2 = variance(data, index, begin, end_, 2);
        if n == 1 || var1 >= var2
            k = 1;
        else
            k = 2;
        end
        index = nth_element(index, begin, begin+floor(n/2), end_, @(a,b)data(a,k)<data(b,k));
        assert(all(data(index(begin:begin+floor(n/2)-1),k) <= data(index(begin+floor(n/2)),k)));
        assert(all(data(index(begin+floor(n/2)+1:end_-1),k) >= data(index(begin+floor(n/2)),k)));
        temp = index(begin+floor(n/2));
        index(begin+floor(n/2)) = index(begin);
        index(begin) = temp;
        nodeslen = nodeslen+1;
        node_index = nodeslen;
        nodes(node_index,1) = k;
        nodes(node_index,2) = data(index(begin),k);
        stacklen = stacklen+1;
        stack(1,stacklen) = begin+floor(n/2)+1;
        stack(2,stacklen) = end_;
        stack(3,stacklen) = node_index;
        stacklen = stacklen+1;
        stack(1,stacklen) = begin+1;
        stack(2,stacklen) = begin+floor(n/2)+1;
        stack(3,stacklen) = node_index;
        if parent ~= 0
            if isnan(nodes(parent,3))
                nodes(parent,3) = node_index;
            else
                nodes(parent,4) = node_index;
            end
        end
    end
    tree.nodes = nodes;
    tree.index = index;
end
function data = filt_by_velocity(data)
    vel = data(:,3);
    upper = mean(vel)+3*std(vel);
    data = data(vel <= upper,:);
end