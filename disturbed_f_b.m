% Lab 3 b) DISTURBED

clear all, clc, clf, close all
format long


function final_main_output = program_b(mass, g, k, y_start, x_start, net_height)

mass = mass  ; 
g = g ; 
k = k ; 
y_start = y_start ; 
x_start  = x_start ; 

%Funktioner

%RK4_Engine
function next_value = RK4_Engine(h, t, value, mass, g,k)

    
    y_bis = @(t, prim_values) (-k*prim_values(2)*sqrt(prim_values(1).^2 + prim_values(2).^2) - mass * g) / mass; 
    x_bis = @(t, prim_values) (-k*prim_values(1)*sqrt(prim_values(1).^2 + prim_values(2).^2)) / mass;
    
    k1 = [value(3); value(4); x_bis(t, [value(3), value(4)]); y_bis(t, [value(3), value(4)])];
    k2 = [value(3) + h/2*k1(3); value(4) + h/2*k1(4); x_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)]); y_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)])];
    k3 = [value(3) + h/2*k2(3); value(4) + h/2*k2(4); x_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)]); y_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)])];
    k4 = [value(3) + h*k3(3); value(4) + h*k3(4); x_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)]); y_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)])];
    
    next_value = value + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end


%Funktion som hittar studspunkterna, via interpolation
%Grad 2
function [first_root,koeff, f_prim] = find_root(guess, x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ;
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ;
    
    %interpolera och definera funktionen
    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ;
    
    % x = -0.5:0.0001:1 ;
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    f_prim = @(x) B(2) + 2*x*B(3);
    first_root = newtons_raphson(guess,f,f_prim); %använd newton för att hitta nollställe!!
    % plot(x, f(x), 'red')
    % grid on
    % hold on
    koeff = [B(1) B(2) B(3)] ;
end


%Funktion som hittar höjden vid x = 0, via interpolation 
%Grad 2
function net_touch = single_interpolation(x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)  
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ;
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ;
    
    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ;
    %x = x_values_:0.0001:x_values_i+1
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    net_touch = f(0);
end


%Funktion som hittar farten i y-led vid studspunkten, via interpolation
%Grad 2
function net_touch = interp_yprim(root,x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)  % annan interpolationsfunktion, till att hitta höjden vid x = 0
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ;
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ;
    
    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ;
    %x = x_values_:0.0001:x_values_i+1
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    net_touch = f(root);
end


%Newtons metod för att hitta nollställen till interpolationspolynomet
function output = newtons_raphson (current_x, f, f_prim)

    xtest = 0 ;
    next_x = 1 ;
    while abs(xtest - current_x) > 10e-8
        next_x = current_x - f(current_x)/f_prim(current_x);
        xtest = current_x;
        current_x = next_x;
    
    end
    output = next_x ;
end

%funktion av programmet i a) that som returnernar höjden över nätet.
    function height = pingpong(speed, mass, g, k, y_start, x_start, net_height)

    % constants 
    t = 0 ;

    h = 0.0002 ; %steplength
  

    %start value
    y_prim_start = 0 ;
    x_prim_start = -speed ;

   
    % intialise values
    u1 = x_start ;
    u2 = y_start ;
    u3 = x_prim_start ;
    u4 = y_prim_start ;

    u = [u1 u2 u3 u4]' ;


    iteration = 8000 ;
    x_values = ones(iteration, 1) ;
    y_values = ones(iteration, 1 ) ;
    y_prim_values = ones(iteration,1);

    root_values = zeros(1,2);
    placement = 1 ;
    test = 0 ;

    for i = 1:iteration
        next_value = RK4_Engine(h, t, u, mass,g,k) ;
        u = next_value ;
        x_values(i) = u(1) ;
        y_values(i) = u(2) ;
        y_prim_values(i) = u(4);

        evaluated_y_value = y_values(i) ;
        evaluated_x_value = x_values(i) ;


        %bounce condition
        if evaluated_y_value < 0

            % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
            [root, koeff, f_prim] = find_root(x_values(i) - 0.1, x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2));

            real_yprim = interp_yprim(root, x_values(i),x_values(i-1),x_values(i-2),y_prim_values(i), y_prim_values(i-1), y_prim_values(i-2));
        


            u(2) = 0 ;
            u(1) = root ;
            u(4) = -real_yprim ;
            x_values(i) = u(1) ;
            y_values(i) = u(2) ;
            root_values(placement) = root;
            placement = placement + 1 ;

        end

        if evaluated_x_value < 0 && test == 0

            height = single_interpolation(x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2));
            

            test = test + 1 ;
        end
    end




end



% Sekantmetoden
    function root = secant(guess_1, guess_2, mass, g, k, y_start, x_start, net_height)

tol = 10e-8;

trial = @(guess) pingpong(guess, mass, g, k, y_start, x_start, net_height) - net_height;

xi = guess_1; % First guess  
xi_minus_1 = guess_2; % Second guess  

f = trial(xi);
g = trial(xi_minus_1);

x1 = 10;
x0 = 0;
A = [];
i = 1;

while abs(guess_1 - guess_2) > tol

    x1 = xi - (f * (xi - xi_minus_1)) / (f - g);
    x0 = xi;
    xi_minus_1 = xi;
    xi = x1;
    f = trial(xi);
    g = trial(xi_minus_1);
    
    A(i) = x1;
    i = i + 1;

    next_value = guess_1 - trial(guess_1) * (guess_1 - guess_2) / (trial(guess_1) - trial(guess_2));
    guess_2 = guess_1;
    guess_1 = next_value;
    root = next_value ;
    
end



end




% second_root = secant(2.8, 4, mass, g, k, y_start, x_start, net_height) ; 
% lola = pingpong(second_root, mass, g, k, y_start, x_start, net_height) ; 


third_root = secant(4, 6, mass, g, k, y_start, x_start, net_height) ;
lola = pingpong (third_root, mass, g, k, y_start, x_start, net_height) ; 
final_main_output = {third_root, lola} ; 
clf, close all 



end 


% ------------------------------- % 

mass = 0.01; 
g = 9.82; 
k = 0.005;

y_start = 0.31; 
x_start = 1.21; 
net_height = 0.119 ; 


dexter = {mass, g, k, y_start, x_start, net_height} ;

iteration = 6 ; 

values = {}; 

for i = 1:iteration


percent = 0.99 ; 
mass = 0.01; 
g = 9.82; 
k = 0.005; 
y_start = 0.31; 
x_start = 1.21; 
net_height = 0.119 ;



dexter = {mass, g, k , y_start, x_start, net_height} ;
texas = dexter{i} * percent ;  %väljer 1 i taget 
dexter{i} = texas ; 


values{i} =     program_b(dexter{:}) ; 
end 

matrix_values = cell2mat(cellfun(@(x) cell2mat(x), values, 'UniformOutput', false))
