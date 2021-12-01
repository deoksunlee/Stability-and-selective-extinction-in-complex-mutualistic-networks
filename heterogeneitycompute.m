function h = heterogeneitycompute(Amatrix,c)
% Atmatrix : NP x NA matrix for mutualistic interaction
% c : competition strength

% h : heterogeneity

[NP NA]=size(Amatrix);
N = [NP NA];
kP = sum(Amatrix,2)';
kA = sum(Amatrix,1);
kmean = [mean(kP) mean(kA)];
k2mean = [mean(kP.*kP) mean(kA.*kA)];
xi = k2mean./(kmean.^2);
ct = (c*N)./(c*N+1-c);
h=sqrt(kmean(1)*kmean(2).*(xi(1)-ct(1))*(xi(2)-ct(2)));
end