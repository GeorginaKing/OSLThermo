function way = randpathAD(origin,target)
v2go = norm(origin-target); % distance between origin and target
vlim = 0.01*v2go; % maximal step is 1% total distance
way = [origin;target];
while max(v2go)>vlim
	vector = way(1:end-1,:)-way(2:end,:);
	v2go = sqrt(vector(:,1).^2 + vector(:,2).^2);
	ix = find(v2go == max(v2go));
	new_point = way(ix,:)-([rand rand].*vector(ix,:));
	way = [way(1:ix,:); new_point; way(ix+1:end,:)];
end

