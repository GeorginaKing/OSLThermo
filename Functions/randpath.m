function out = randpath(origin, target);
% written by Benny Guralnik (benny.guralnik@gmail.com)
% Usage:
% returns a random monotonous path from [origin] to [target] using steps
% no bigger than [v_min]
% Input arguments:
% origin and target are 2D vectors of the form [x y]
% Output arguments:
% out is a vector containting cosecutive [x y] pair in each row, i.e.
% [x1 y1; x2 y2; ... ; xn yn;]
% Example:
% a = randpath([0 0],[1 1])
% plot(a(:,1),a(:,2))
v_min = 0.0001; % minimal displacement (out of 1)
path = []; % path history
current_position = origin;
v_togo = target-origin;
while norm(v_togo) > v_min;
v = [target - current_position].*[rand rand]; %random displacement
if norm(v) > v_min % break a large displacement further down by recursion
path=[path; randpath(current_position, current_position + v)];
end
path=[path; current_position + v ]; % save current position
v_togo = v_togo - v;
current_position = current_position + v;
end
out=path;