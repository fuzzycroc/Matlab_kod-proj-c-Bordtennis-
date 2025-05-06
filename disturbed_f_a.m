% Lab 3 a) STÖRNING 

%main_code

clear all, clc, clf, close all
format long


function final_output = program_a(mass, g, k_x, k_y,  velocity, y_start,x_start)



% constants 
length = 1.21 ; 
height = 0.119 ;  
t = 0 ; 
mass = mass  ; 

y_start = y_start ; 
x_start = x_start ; 


h = 0.0002 ;



y_prim_start = 0 ; 
x_prim_start = -velocity ; 
iteration = 8000 ;


%functions 

%RK4_Engine
function next_value = RK4_Engine(h, t, value, mass, g, k_x, k_y)

    
    
    y_bis = @(t, prim_values) (-k_y*prim_values(2)*sqrt(prim_values(1).^2 + prim_values(2).^2) - mass * g) / mass; 
    x_bis = @(t, prim_values) (-k_x*prim_values(1)*sqrt(prim_values(1).^2 + prim_values(2).^2)) / mass;
    
    k1 = [value(3); value(4); x_bis(t, [value(3), value(4)]); y_bis(t, [value(3), value(4)])];
    k2 = [value(3) + h/2*k1(3); value(4) + h/2*k1(4); x_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)]); y_bis(t + h/2, [value(3) + h/2*k1(3), value(4) + h/2*k1(4)])];
    k3 = [value(3) + h/2*k2(3); value(4) + h/2*k2(4); x_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)]); y_bis(t + h/2, [value(3) + h/2*k2(3), value(4) + h/2*k2(4)])];
    k4 = [value(3) + h*k3(3); value(4) + h*k3(4); x_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)]); y_bis(t + h, [value(3) + h*k3(3), value(4) + h*k3(4)])];
    
    next_value = value + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end

%Interpolerar mellan x-och y-värdena vid studs, avsedd för att hitta
%studspunkten
%grad 2
function [first_root,koeff, f_prim] = find_root(guess, x_values_i, x_values_i_1, x_values_i_2, y_values_i, y_values_i_1, y_values_i_2)
    interpolate_x_values = [x_values_i x_values_i_1 x_values_i_2]' ; 
    interpolate_y_values = [y_values_i y_values_i_1 y_values_i_2]' ; 
    
    %interpolera och definera funktionen
    A = ones(3,3);
    A(:, 2) = interpolate_x_values ;
    A(:, 3) = interpolate_x_values.^2 ;
    
    B = A\interpolate_y_values ; 
    
    x = -0.5:0.0001:1 ;
    f = @(x) B(1) + B(2)*x + B(3)*x.^2 ;
    f_prim = @(x) B(2) + 2*x*B(3);
    first_root = newtons_raphson(guess,f,f_prim); %använd newton för att hitta nollställe!!
    % plot(x, f(x), 'red')
    % grid on
    % hold on
    koeff = [B(1) B(2) B(3)] ; 
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


%Newtons metod för att hitta nollställen till interpolationspolynomet
function output = newtons_raphson(current_x, f, f_prim)

p = 8; 

x1 = current_x;  
f_val = f(x1);
fprim_val = f_prim(x1);

xi = 10; 
xtest = 0;
i = 1;

while abs((xi - xtest) / xi) > 10^-p   
    
    xi = x1 - f_val / fprim_val;
    xtest = x1;
    x1 = xi;
    f_val = f(x1); 
    fprim_val = f_prim(x1);
    i = i + 1;
end
    
output = x1 ; % Final computed root



end


% % Interpolationsfunktion för att hitta felet när vi interpolerar mellan x-
% % och y-värden vid studs
% %grad 1 
% function root = interp_fel(guess, x_values_i, x_values_i_1,y_values_i, y_values_i_1) 
%     interpolate_x_values = [x_values_i x_values_i_1]';
%     interpolate_y_values = [y_values_i y_values_i_1]';
%     A = ones(2,2) ;
%     A(:,2) = interpolate_x_values ;
%     B = A\interpolate_y_values;
%     f = @(x) B(1) + B(2)*x;
%     f_prim = @(x) B(2);
%     root = newtons_raphson(guess,f,f_prim);
% 
% 
% end

    % intialise values 
    u1 = x_start ;
    u2 = y_start ;
    u3 = x_prim_start ;
    u4 = y_prim_start ;
    
    t = 0; 
    u = [u1 u2 u3 u4]' ;
    %iteration = 8000 
    x_values = ones(iteration, 1) ;
    y_values = ones(iteration, 1 ) ;
    y_prim_values = ones(iteration,1) ; 

  

    root_values = zeros(1,2);
    placement = 1 ;
    x_counter = 1;

    for i = 1:iteration
        next_value = RK4_Engine(h,t,u, mass, g, k_x, k_y) ;
        u = next_value ;
        x_values(i) = u(1) ;
        y_values(i) = u(2) ;
        y_prim_values(i) = u(4);
        evaluated_y_value = y_values(i) ;
        evaluated_x_values = x_values(i);


        %bounce condition
        if evaluated_y_value < 0 % Hitta nollställe och ändra lite på vektorn utifrån nollstället

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
        if evaluated_x_values < 0 && x_counter == 1
            net_touch1 = single_interpolation(x_values(i),x_values(i-1),x_values(i-2),y_values(i), y_values(i-1), y_values(i-2));
           
            x_counter = x_counter + 1;
        end
 
    end



final_output = {root_values(1), root_values(2), net_touch1} ;





end 







% mass = 0.01; 
% g = 9.82 ; 
% k_x = 0.005 ;
% k_y = 0.005 ;
% velocity = 4; 
% y_start = 0.31; 
% x_start = 1.21 * 0.99; 


% dexter = {mass, g, k_x, k_y ,velocity, y_start, x_start} ;


% program_a(mass, g, k_x, k_y ,velocity, y_start, x_start) ; 

%orginala värde
%  {[0.351637568704311]}    {[-0.430633009491603]}    {[0.183655472911037]}


% ------------------------- %
% mass * 1.01 
% {[0.350396869361334]}    {[-0.437249940973851]}    {[0.183750766770818]}

% g * 1.01 
%  {[0.355308288276389]}    {[-0.425276604936141]}    {[0.184639002259837]}

% k_x * 1.01 
%  {[0.353408208136643]}    {[-0.424937519505446]}    {[0.184279371316807]}

% k_y * 1.01 
% {[0.351113438486938]}    {[-0.429699097534503]}    {[0.182910914859762]}

% veocity * 1.01 
%{[0.344256566493537]} {[-0.441288336212319]} {[0.181568363568636]}

% y_start * 1.01 
% {[0.348012900370209]}    {[-0.434865435953440]}    {[0.184457665513697]}


% x_start * 1.01 
%    {[0.363737568704311]}    {[-0.418533009491651]}    {[0.185408837731868]}


%------------------------------%
% mass * 0.99
% {[0.352899553131431]}    {[-0.423950639670327]}    {[0.183531964259155]}


% g * 0.99
%  {[0.347916578467645]}    {[-0.436075208957321]}    {[0.182614986892128]}


% k_x * 0.99
% {[0.349857779783480]}    {[-0.436367971588925]}    {[0.183006832492286]}

% k_y * 0.99
% {[0.352161336709142]}    {[-0.431567739646093]}    {[0.184401757249220]}

% veocity * 0.99
%{[0.359039240736829]}    {[-0.419702114465735]}    {[0.185607175553185]}

% y_start * 0.99
%    {[0.355286447088528]}    {[-0.426141006473486]}    {[0.182830737117699]}


% x_start * 0.99
%    {[0.339537568704311]}    {[-0.442733009491584]}    {[0.181560742285943]}


% Define origin and perturbation vectors
origin = [0.351637568704311, -0.430633009491603, 0.183655472911037];

% Define all +1% and -1% variations as rows in matrices
plus_variations = [
    0.350396869361334, -0.437249940973851, 0.183750766770818;  % mass1
    0.355308288276389, -0.425276604936141, 0.184639002259837;  % g1
    0.353408208136643, -0.424937519505446, 0.184279371316807;  % kx1
    0.351113438486938, -0.429699097534503, 0.182910914859762;  % ky1
    0.344256566493537, -0.441288336212319, 0.181568363568636;  % velocity1
    0.348012900370209, -0.434865435953440, 0.184457665513697;  % ystart1
    0.363737568704311, -0.418533009491651, 0.185408837731868   % xstart1
];

minus_variations = [
    0.352899553131431, -0.423950639670327, 0.183531964259155;  % mass9
    0.347916578467645, -0.436075208957321, 0.182614986892128;  % g9
    0.349857779783480, -0.436367971588925, 0.183006832492286;  % kx9
    0.352161336709142, -0.431567739646093, 0.184401757249220;  % ky9
    0.359039240736829, -0.419702114465735, 0.185607175553185;  % velocity9
    0.355286447088528, -0.426141006473486, 0.182830737117699;  % ystart9
    0.339537568704311, -0.442733009491584, 0.181560742285943   % xstart9
];

param_names = {'mass', 'g', 'k_x', 'k_y', 'velocity', 'y_start', 'x_start'};

% Loop through each parameter and compare deviations

iteration = 1 ; 



lord_of_rings =  [] ; 



for i = 1:length(param_names)
    value1 = abs(origin(iteration) - plus_variations(i,iteration));   % abs difference for +1%
    value2 = abs(origin(iteration) - minus_variations(i,iteration));  % abs difference for -1%
    
   
    if value1 > value2
        disp(['For ', param_names{i}, ', value1 (+1%) caused more change.']);
        lord_of_rings(i) = value1 ; 
    else
        disp(['For ', param_names{i}, ', value2 (-1%) caused more change.']);
        lord_of_rings(i) = value2 ; 
    end
end

sum(lord_of_rings) ; 
lord_of_rings
%Första_studs =  0.030437444218725
% Andra_studs = 0.046317159582772 
% Nät_höjd  =   0.007565495188558

