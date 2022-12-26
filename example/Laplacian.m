function L = Laplacian(file)
graph = readmatrix(file);
nv = graph(1,1);
ne = graph(1,2);
L = zeros(nv,nv);
for i = 2:ne+1
    L(graph(i,1),graph(i,2)) = -graph(i,3);
    L(graph(i,2),graph(i,1)) = L(graph(i,1),graph(i,2));
    L(graph(i,1),graph(i,1)) = L(graph(i,1),graph(i,1)) + graph(i,3);
    L(graph(i,2),graph(i,2)) = L(graph(i,2),graph(i,2)) + graph(i,3);
end
end
