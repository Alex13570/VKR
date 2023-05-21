function [P] = applied_force(t, n,cT)

% if t*4 < 1
%     k = t*4;
% else 
%     k = 1;
% end

P = [
    zeros(n,1);
    ones(n,1)*cos(cT*t);
    zeros(n,1) ;
    zeros(n,1);
];

end