function Mtsp = convrtTspCelltoMatrix(Mtspcell);

nrpts = length(Mtspcell);
maxnsp = max(celleval('length', Mtspcell));

Mtsp = zeros(maxnsp,nrpts);
for j = 1:nrpts
    tsp = Mtspcell{j};
    Mtsp(1:length(tsp),j) = tsp;
end
