% Kontroll noggrannhetsordning

clear all, clc, clf, close all
format long

% constats 
mass = 0.01 ; 
k = 0.005 ;
length = 1.21 ; 
height = 0.119 ; 
g = 9.82 ;
t = 0 ;

h = 0.002 ; %steplength 

%start values 

y_start = 0.31 ; 
x_start = 1.21 ; 
y_prim_start = 0 ; 
x_prim_start = -4 ; 

%functions 

%RK4_Engine

function next_value = RK4_Engine(h, t, value)

    mass = 0.01; 
    k = 0.005;
    g = 9.82;
    
    y_bis = @(t, prim_values) (-k*prim_values(2)*sqrt(prim_values(1).^2 + prim_values(2).^2) - mass * g) / mass; 
    x_bis = @(t, prim_values) (-k*prim_values(1)*sqrt(prim_values(1).^2 + prim_values(2).^2)) / mass;
    
    k1 = [value(3); value(4); x_bis(t, [value(3), value(4)]); y_bis(t, [value(3), value(4)])];
    k2 = [value(3) + h/2*k1(3); value(4) + h/2*k1(4); x_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)]); y_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)])];
    k3 = [value(3) + h/2*k2(3); value(4) + h/2*k2(4); x_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)]); y_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)])];
    k4 = [value(3) + h*k3(3); value(4) + h*k3(4); x_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)]); y_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)])];
    
    next_value = value + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end

u1 = x_start ;
u2 = y_start ;
u3 = x_prim_start ; 
u4 = y_prim_start ; 

u = [u1 u2 u3 u4]' ; 


% next_value = u + h*RK4_Engine(h,t,u) ;


iteration = 40 ; 
x_values = ones(iteration, 1) ;
y_values = ones(iteration, 1 ) ; 



% % change bounce condition 
root = 0 ; 
new_y_prim = 0 ; 
root_values = zeros(1,2); 
placement = 1 ;
x_counter = 1;
for i = 1:iteration


    next_value = RK4_Engine(h,t,u) ;
    u = next_value ;
    x_values(i) = u(1) ;
    y_values(i) = u(2) ;

    evaluated_y_value = y_values(i) ;
    evaluated_x_values = x_values(i);
end
%För att jämföra ändra steglängden enligt nedan tänk på att ändra
%motsvarande iterationer. 
slut_varde = x_values(end) 

% h = 0.001, x = 0.913230759127980
% h = 0.002, x = 0.913230759131567
% h = 0.0005, x = 0.913230759127757

fel1 = abs(0.913230759127980 - 0.913230759131567)
fel2 = abs(0.913230759127980 - 0.913230759127757)
nog = fel1/16 % vi ser att nog = fel2