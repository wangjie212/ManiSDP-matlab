function problem = sedumi2copt(A, b, c, K)
% Convert a SeDuMi-format problem into a COPT model struct
    
    ind = 1:K.f;
    for k = 1:length(K.s)
        a = sum(K.s(1:k-1).^2) + K.f;
        for i = 1:K.s(k)
            ind = [ind (i:K.s(k)) + (i-1)*K.s(k) + a];
        end
    end

    model.A = A(:, ind);
    model.b = b;
    model.c = c(ind);
    model.K = K;
    model.objsen = 'maximize';
    problem.conedata = model;

end
