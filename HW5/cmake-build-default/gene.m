data = textread('bestGene.txt', '', 'delimiter', ' ');
data = data(1600,1:484);
A = "";
for i = 1:484
    A = strcat(A, ",", num2str(data(i)));
end
A