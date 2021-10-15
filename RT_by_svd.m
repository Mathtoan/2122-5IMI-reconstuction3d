function [R, T] = RT_by_svd(p, q)
    P = zeros(size(p));
    Q = zeros(size(q));


    for i=1:length(p)
        P(i,:) = p(i,:) - mean(p);
        Q(i,:) = q(i,:) - mean(q);
    end

    P = transpose(P);
    Q = transpose(Q);

    % Calcul M et de sa SVD
    M = P*transpose(Q);
    [U,~,V] = svd(M);

    % Matrice de rotation
    Id = eye(3);
    Id(end:end) = det(V*transpose(U));
    R = V*Id*transpose(U);

    T = transpose(mean(q)) - R*transpose(mean(p));
end