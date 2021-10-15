clear variables;
close all;

% Creation de P et Q
tic
p_grossier = [-208.710999,  42.494999, 126.382004;
              -248.830994, 190.302994, 255.768997;
              -339.148010, 190.569000, 103.722000;
              -158.994003,  12.918700, 184.602005];

q_grossier = [-132.136993,  36.402302, 178.623001;
              - 43.392200, 182.332993, 288.166992;
              -220.121994, 183.210999, 279.062012;
              - 56.704700,   6.993410, 169.884995];
 
[R_grossier, T_grossier] = RT_by_svd(p_grossier, q_grossier);
toc

% Recalage grossier
tic
set0 = load('set0.xyz');
set1 = load('set1.xyz');

set0 = transpose(set0(1:end, 1:3));
set1 = transpose(set1(1:end, 1:3));


[col,ligne] = size(set0);
set0_grossier = zeros(col, ligne);

for i=1:ligne
    set0_grossier(:,i) = R_grossier*set0(:,i) + T_grossier;
end
toc
temp = set0_grossier';
save 'set0_grossier.xyz' temp -ascii
% ICP

it = 10;
delta = 0.5;
tic
dt = DelaunayTri(transpose(set1));
toc

tic
%% 
set0_int = set0_grossier;
for k=1:it
    [pid, d] = nearestNeighbor(dt, transpose(set0_int));
    indice = find(d<= delta);
    
    p_int = transpose(set0_int(:, indice));
    q_int = transpose(set1(:, pid(indice)));
%     for i=1:ligne
%         if d(i) <= delta
%             p_int = [p_int; transpose(set0_int(:, i))];
%             q_int = [q_int; transpose(set1(:, pid(i)))];
%             c = c + 1;
%         end
%     end
    disp(length(indice));
    [R_int, T_int] = RT_by_svd(p_int, q_int);
    
    for j=1:ligne
        set0_int(:,j) = R_int*set0_int(:,j) + T_int;
    end
end
toc

set0_int = set0_int';

save 'set0_int_delta05_10.xyz' set0_int -ascii