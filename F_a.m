% give in a

% Lab 3 a) 

clear all, clc, clf, close all
format long

% Konstanter 
mass = 0.01 ; 
k = 0.005 ;
length = 1.21 ; 
height = 0.119 ; 
g = 9.82 ; 
t = 0 ;

% steglängd

h = 0.0002 ;

%Startvärden
y_start = 0.31 ; 
x_start = 1.21 ; 
y_prim_start = 0 ; 
x_prim_start = -4 ; % Passar det definierade koordinatsystemet
iteration = 8000 ; % Kan väljas fritt, ju fler iterationer desto längre kör programmet


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

%Interpolerar mellan x-och y-värdena vid studs, avsedd för att hitta
%studspunkten
%grad 2
function root = find_root(guess, x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ; 
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ; 
    
    %Ställ upp matrisekvationen
    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ; 
    
    x = -0.5:0.0001:1 ;
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    f_prim = @(x) B(2) + 2*x*B(3);
    root = newtons_raphson(guess,f,f_prim); %använd newton för att hitta nollställe!!

end

% Interpolationsfunktion för att hitta felet när vi interpolerar mellan x- och y-värden vid studs
%grad 1 
function root = interp_fel(guess, x_values_i, x_values_i_1,y_values_i, y_values_i_1) 
    interpolate_x_values = [x_values_i x_values_i_1]';
    interpolate_y_values = [y_values_i y_values_i_1]';
    A = ones(2,2) ;
    A(:,2) = interpolate_x_values ;
    B = A\interpolate_y_values;
    f = @(x) B(1) + B(2)*x;
    f_prim = @(x) B(2);
    root = newtons_raphson(guess,f,f_prim);

    
end

% annan interpolationsfunktion, till att hitta höjden vid x = 0
% grad 2
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

%interpolationsfunktion för att hitta felet när vi interpolerar vid x = 0. sänkt gradtal. 
%grad 1
function net_touch = net_interp_error(x_values_i, x_values_i_1, y_values_i, y_values_i_1) 
    interpolate_x_values = [x_values_i x_values_i_1]';
    interpolate_y_values = [y_values_i y_values_i_1]';
    A = ones(2,2);
    A(:, 2) = interpolate_x_values ;
    B = A\interpolate_y_values;
    f = @(x) B(1) + B(2)*x;
    net_touch = f(0);
end

%interpolerar mellan yprim- och x-värdena värden vid studserna, används för att hitta hastighet vid studs
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

%interpolationsfunktion för att hitta felet när vi interpolerar mellan yprim- och x-värdena vid studsarna. sänkt gradtal.
% interpolerar vid studs
function bounce_touch = single_interp_yprim(root,x_values_i, x_values_i_1, y_values_i, y_values_i_1)  
    interpolate_x_values = [x_values_i x_values_i_1]';
    interpolate_y_values = [y_values_i y_values_i_1]';
    A = ones(2,2);
    A(:, 2) = interpolate_x_values ;
    B = A\interpolate_y_values;
    f = @(x) B(1) + B(2)*x;
    bounce_touch = f(root);
end

%Newtons metod för att hitta nollställen till interpolationspolynomet
function output = newtons_raphson(current_x, f, f_prim)

p = 8; 

x1 = current_x;  
f_val = f(x1);
fprim_val = f_prim(x1);

xi = 10; 
xtest = 0;
i = 1;
A = [];

while abs(xi - xtest) > 10^-p   
    A(i) = x1;
    xi = x1 - f_val / fprim_val;
    xtest = x1;
    x1 = xi;
    f_val = f(x1); 
    fprim_val = f_prim(x1);
    i = i + 1;
end
    
output = x1 ; % Final computed root

a = length(A);
i = 1;
B = [];

for c = 0:a - 2
   B(i) = abs(A(c+1) - A(c+2)); % Calculates error terms for each iteration, places them in a vector
   i = i + 1;
end

P = 2; % test P value 


for j = 1:length(B) - 1
    C = B(j + 1) / (B(j) ^ P) ;% Checks convergence constant
end

end


    % initialvärden 
    u1 = x_start ;
    u2 = y_start ;
    u3 = x_prim_start ;
    u4 = y_prim_start ;
    
    t = 0; 
    u = [u1 u2 u3 u4]' ;
    next_value = RK4_Engine(h,t,u) ; % Börjar räkna på DE:n
    x_values = ones(iteration, 1) ;
    y_values = ones(iteration, 1 ) ;
    y_prim_values = ones(iteration,1) ; 

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
        y_prim_values(i) = u(4);
        evaluated_y_value = y_values(i) ;
        evaluated_x_values = x_values(i);


        %bounce condition
        if evaluated_y_value < 0 % Hitta nollställe och ändra lite på vektorn utifrån nollstället
            
             % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
            root = find_root(x_values(i) - 0.1, x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2));
            %root = interp_fel(x_values(i) - 0.1,x_values(i),x_values(i-1),y_values(i), y_values(i-1));
           
            %Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
            real_yprim = interp_yprim(root, x_values(i),x_values(i-1),x_values(i-2),y_prim_values(i), y_prim_values(i-1), y_prim_values(i-2)) ;
            %real_yprim = single_interp_yprim(root, x_values(i),x_values(i-1),y_prim_values(i), y_prim_values(i-1)) ; 
            

            u(2) = 0 ;
            u(1) = root ;
            u(4) = -real_yprim ;
            x_values(i) = u(1) ;
            y_values(i) = u(2) ;

            root_values(placement) = root;
            placement = placement + 1 ;

        end
        if evaluated_x_values < 0 && x_counter == 1

            % Kommentera bort det interpolationspolynom som inte ska användas, för felet, ta skillnaden.
            net_touch = single_interpolation(x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2));
            %net_touch = net_interp_error(x_values(i),x_values(i-1),y_values(i), y_values(i-1));

            x_counter = x_counter + 1;
        end
    end


% Plottar bollbanan

% plot(x_values, y_values, color = 'blue') ;
%     grid on
%     hold on
% 
% % Plottar nätet och bordet
% board_length = 1.21 ; 
% net_height = 0.119 ; 
% net_values = 0:0.001:net_height ;
% board_values = -board_length:0.001:board_length ; 
% net = @(net_values) 0.*net_values  ; 
% board = @(board_values) 0.*board_values  ; 
% plot(net(net_values), net_values, 'red')
% hold on 
% plot(board_values, board(board_values), 'red')
% hold on

%Svaret vad gäller studspunkterna
root_values(1)
root_values(2)
%Svaret vad gäller höjden vid nätet
net_touch  
