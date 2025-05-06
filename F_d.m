% give in d

% del d)

clear all, clc, clf, close all
format long

% steglängden/tiden per iteration
h = 0.0002;
    

%Funktioner
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

%Interpolerar mellan y- och x-värdena vid studsarna
% grad 2
function root = find_root(guess, x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)

    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ; 
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ; 
    
 
    A = ones(3,3) ;
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ; 
    
    %x = -0.5:0.0001:1 ;
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    f_prim = @(x) B(2) + 2*x*B(3);
    root = newtons_raphson(guess,f,f_prim) ;
end


%interpolerar mellan yprim- och x-värdena värden vid studserna, används för att hitta hastighet vid
%studs, Används även till interpolation mellan tidsvärdena och x-värdena
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


%Newtons metod för att hitta nollställen till vissa av interpolationspolynomen
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

% Skickar ut index för antal iteerationer vid varje studs och tiden vid
% varje studs
function [B,a] = find_height(angle, interpolations_metod)

% constants 
mass = 0.01 ; 
k = 0.005 ;
length = 1.21 ; 
height = 0.119 ; 
g = 9.82 ;
total_time = 0 ; 
speed_start = -10 ; 

jumanji = 0 ; 

h = 0.0002; % steglängden/tiden per iteration


%Starvärden
y_start = 0.31 ; 
x_start = 1.21 ; 
y_prim_start =  sin(angle)*speed_start; 
x_prim_start = cos(angle)*speed_start; 


% intialise values 
u1 = x_start ;
u2 = y_start ;
u3 = x_prim_start ; 
u4 = y_prim_start ; 
u = [u1 u2 u3 u4]' ; 

x_values = [] ;
y_values = [] ; 
y_prim_values= [];

root_values = zeros(1,2); 
placement = 1 ;
B = [] ; 

i = 1 ; 
j = 1 ;
a = zeros(1,2);


while placement <= 2

    next_value = RK4_Engine(h,total_time,u) ;
    u = next_value ;
    x_values(i) = u(1) ;
    y_values(i) = u(2) ;
    y_prim_values(i) = u(4);

    evaluated_y_value = y_values(i) ;
    evaluated_x_value = x_values(i) ;

    %bounce condition 
    if evaluated_y_value < 0 % Hitta nollställe och ändra lite på vektorn utifrån nollstället
  
        % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
        root = find_root(x_values(i) - 0.1, x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2)) ;
      
        % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
        real_yprim = interp_yprim(root,x_values(i),x_values(i-1),x_values(i-2),y_prim_values(i), y_prim_values(i-1), y_prim_values(i-2));
       

        % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
        root_time = interp_yprim(root,x_values(i),x_values(i-1),x_values(i-2),h*(i),h*(i-1),h*(i-2));
    
        u(2) = 0 ;
        u(1) = root ;
        u(4) = -real_yprim;
        x_values(i) = u(1) ;
        y_values(i) = u(2) ; 

        a(j) = i;
        B(j) = root_time;
        root_values(placement) = root  ;
        placement = placement + 1 ; 
        j = j+1;

     end 
       
     i = i + 1 ;
    
end 

%Plottar
% net_values = -1.21:0.01:0 ;
% test_values = zeros(1,numel(net_values)) ;
% test_values(end) = 0.119 ;
% plot(net_values, test_values)
% hold on
% plot(x_values, y_values, color = 'blue') ; 
% grid on

end 

[root_time,Index] = find_height(0.853126126231213);

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
% xlabel('x-position')
% ylabel('y-position')



% svar:
real_time = h*Index(2) - ((h*Index(2) - root_time(2))+(h*Index(1) - root_time(1))) %  0.935341533875251

