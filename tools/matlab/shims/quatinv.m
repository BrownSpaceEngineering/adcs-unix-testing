function r = quatinv(q)
% Shim for the Aerospace Toolbox function: conjugate / |q|^2 (1x4, scalar first).
r = [q(1), -q(2), -q(3), -q(4)] / sum(q.^2);
end
