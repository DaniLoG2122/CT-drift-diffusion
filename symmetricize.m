function y=symmetricize(sol)

% Function to symmetricize solution

xp = length(sol.x);

symsol = sol;
symsol.sol(:, 1:xp, :) = sol.sol(:,:,:);

symsol.x(1:xp) = sol.x;
for ii=2:xp
symsol.sol(:, xp+ii-1, :) = symsol.sol(:, xp-ii+1, :);
symsol.x(ii+xp-1) = 2*sol.x(end) - sol.x(xp-ii+1);
end
symsol.n(end, 1:xp) = sol.n;
symsol.n(end, xp+1:2*xp) = fliplr(sol.n(end, :));

symsol.p(end, 1:xp) = sol.p;
symsol.p(end, xp+1:2*xp) = fliplr(sol.p(end, :));


symsol.V(end, 1:xp) = sol.V;
symsol.V(end, xp+1:2*xp) = fliplr(sol.V(end, :));
y=symsol;
assignin('base', 'symsol', symsol);

end