function [y] = next_step(t, A, Ux, Ux_1, mas, cT)
n = length(Ux)/4;
y = (mas-A)\(+2*mas*Ux - mas*Ux_1-[applied_force(t, n, cT)]);
end