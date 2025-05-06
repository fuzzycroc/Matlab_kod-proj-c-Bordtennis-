% give in c
% Lab 3 c) 

clear all, clc, clf, close all
format long


function end_height = find_height(angle)

% constants 
mass = 0.01 ; 
k = 0.005 ;
g = 9.82 ;
t = 0 ; 
speed_start = -10 ; 


h = 0.00004 ; 

%start values 
y_start = 0.31 ; 
x_start = 1.21 ; 
y_prim_start =  sin(angle)*speed_start; 
x_prim_start = cos(angle)*speed_start; 

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

%functin that finds the roots by interpolating values 
function [first_root,koeff, f_prim] = find_root(guess, x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ; 
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ; 
    
 
    A = ones(3,3) ;
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ; 
    
   %x = -0.5:0.0001:1 ;
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    f_prim = @(x) B(2) + 2*x*B(3);
    first_root = newtons_raphson(guess,f,f_prim) ;
    koeff = [B(1) B(2) B(3)] ; 
end


%Newtons metod för att hitta nollställen till interpolationspolynomet
function output = newtons_raphson (current_x, f, f_prim )
xtest = 0 ; 
next_x = 1 ; 
    while abs(xtest - current_x) > 10e-8
        next_x = current_x - f(current_x)/f_prim(current_x);
        xtest = current_x;
        current_x = next_x;
    
    end
output = next_x ; 
end



%interpolerar mellan yprim- och x-värdena värden vid studserna, används för att hitta hastighet vid
%studs
% grad 2 
function bounce_touch = interp_yprim(root,x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2) 
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ; 
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ;

    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ;
    %x = x_values_:0.0001:x_values_i+1 
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    bounce_touch = f(root);
end





% intialise values 
u1 = x_start ;
u2 = y_start ;
u3 = x_prim_start ; 
u4 = y_prim_start ; 
u = [u1 u2 u3 u4]' ; 



x_values = [] ;
y_values = [] ; 
y_prim_values = [];


root_values = zeros(1,2); 
placement = 1 ;

i = 1 ; 


while placement <= 2

    next_value = RK4_Engine(h,t,u) ;
    u = next_value ;
    x_values(i) = u(1) ;
    y_values(i) = u(2) ;
    y_prim_values(i) = u(4);

    evaluated_y_value = y_values(i) ;
   

    %bounce condition 
    if evaluated_y_value < 0 % Hitta nollställe och ändra lite på vektorn utifrån nollstället, välj någon av interpolationsfunktionerna, skillnaden på värdet ger felet. 
        [root, koeff, f_prim] = find_root(x_values(i) - 0.1, x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2)) ;
       

        real_yprim = interp_yprim(root, x_values(i),x_values(i-1),x_values(i-2),y_prim_values(i), y_prim_values(i-1), y_prim_values(i-2)) ;
    

        u(2) = 0 ; 
        u(1) = root ; 
        u(4) = -real_yprim ; 
        x_values(i) = u(1) ;
        y_values(i) = u(2) ; 

        root_values(placement) = root  ;
        placement = placement + 1 ; 


     end 
       
     i = i + 1 ;
    
   


end 


end_height = root_values(2) ;



end 



%------------------%

function root = secant(guess_1, guess_2)

tol = 10e-8 ; 


trial = @(guess) find_height(guess) + 1.21 ; 

while abs(guess_1 - guess_2) > tol 

next_value = guess_1 -  trial(guess_1) * (guess_1 - guess_2) / (trial(guess_1)-trial(guess_2)) ;

guess_2 = guess_1;
guess_1 = next_value;

root = next_value ; 

end 
end 



% board_length = 1.21 ; 
% net_height = 0.119 ; 
% %plots net
% net_values = 0:0.001:net_height ;
% board_values = -board_length:0.001:board_length ; 
% net = @(net_values) 0.*net_values  ; 
% board = @(board_values) 0.*board_values  ; 
% plot(net(net_values), net_values, 'red')
% hold on 
% plot(board_values, board(board_values), 'red')
% hold on
% grid on 


%--------------------------% 
%hittar vilket vinkel fungerar 
clf,close all

% v = -1:0.05:2;               
% nomad = @(input) find_height(input) + 1.21;  
% 
% nomad_values = zeros(size(v));  
% 
% for i = 1:length(v)
%     nomad_values(i) = nomad(v(i));  
% end
% 
% plot(v, nomad_values, 'r');   
% grid on;                      
% xlabel('radianer (rad) ')
% ylabel('slut-bord position')

% svar: hittar vilken vinkel fungerar
vinkel = secant(pi/3,pi/4) %    0.852908753924230
