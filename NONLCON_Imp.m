function [C,Ceq,DC,DCeqt] = NONLCON_Imp(x,e,a,t,w,f,c,ACoef,UCoef,cCoef)

C = [];

e(1) = x(1);
a(1) = x(2);
t(1) = x(3);
%p(1) = x(4);
%p(2) = x(5);

% Vehicle Attributes
Ceq(1) = 1000 ./ ( e(1) - ACoef(1) ) ...
    - ( ACoef(2) ...
    + ACoef(3) .* exp( -a(1) ) ...
    + ACoef(4) .* w(1)...
    + ACoef(5) .* w(1).* a(1) ...
    + ACoef(6) .* t(1) ...
    + ACoef(7) .* t(1) .* a(1)^2 );

if nargout > 2
    
    DC = [];
    
    % Jacobian
    DCeq = zeros(1,5);
    
    % Design constraint
    DCeq(1,1) = -1000/((e(1) - ACoef(1))^2);
    DCeq(1,2) = ACoef(3)*exp(-a(1)) ...
        - ACoef(5)*w(1) ...
        - 2*ACoef(7)*a(1)*t(1);
    DCeq(1,3) = -ACoef(6) - ACoef(7)*a(1)^2;
    
    DCeqt = DCeq';
end

% global Safety;
% 
% if Safety == 1
%     
%     Safety = 0;
%     
%     CheckGrad(@(x)NONLCONZ(x,e,a,t,w,f,c,ACoef,UCoef,cCoef), x, 3);
%     
%     Safety = 1;
%     
% end

end
