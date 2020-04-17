%%
clear
clc

%Rotor
brake_radius_min = 77.7/1000; %m
brake_radius_max = 330/1000; %m
brake_thcknss_min = 20/1000; %m
brake_thcknss_max = 30/1000; %m

nmbr_pads = 8;
torque_needed = 5279.21; %Nm (total)
torque_per_pad = torque_needed / nmbr_pads; %Nm

I = 196.02 / 4; %kg m^2

%%
%Material Properties
c = 490; %specific heat J/kg-K 
h = 41; %heat transfer coef watts/m^2-K 
density = 7500; %kg/m^3
thrml_cndctvty = 41; %W/m-K
thrml_dffsvty = 11; %m^2/s


Ti = 23 + 273; %room temp in Kelvin
v = 160 / 3.6; %m/s


%%
%Pads
P_max = 6*10^6; %pa
mew = .5; 
theta_min = 20; %degrees
theta_max = 90; %degrees
T_max = 800 + 273; %Kelvin
pad_radius_min = 80/1000; %m
pad_density = 1400; %kg/m^3

%%
%calculations
R0 = brake_radius_min; 
c=1;
n=1;
u=1;
while R0 <= brake_radius_max
    
    thcknss = brake_thcknss_min;
    while thcknss <= brake_thcknss_max
        
        volume = pi * (R0^2) * thcknss;
        SurfaceArea = pi * (R0^2);
        m = density * volume;
        w = v/R0; %angular velocity in radians/second 
        
        Ri = pad_radius_min; 
        while Ri <= R0
            
            theta  = theta_min;
            while theta <= theta_max
     
                theta_radians = theta * (pi/180); %radians
                area_pad = (.5*(R0^2)*theta_radians) - (.5*(Ri^2)*theta_radians);
            
                area_touched = 2*(.5*(R0^2)*2*pi) - (.5*(Ri^2)*2*pi);
                e = exp(-h*area_touched/(c*m)*150);
                T = Ti + ((I*(w^2))/(2*c*m*(1-e))); 
                
                Temp_pad = -.5 * theta_radians * mew * P_max * Ri * ( Ri^2 - R0^2);
                
                braking_torque = .5 * area_pad * P_max * ((R0+Ri)/2);
                
                if Temp_pad <= T_max
                
                if braking_torque >= torque_per_pad
                    
                    if T <= T_max * 1.05
            
                        k = area_pad/SurfaceArea;
                        T_final = T * exp(-k * 150); 
            
                        if T_final <= T_max * 1.05
                
                            Rotor(1,n) = R0;
                            Rotor(2,n) = thcknss;
                            Rotor(3,n) = Ri;
                            Rotor(4,n) = theta;
                            Rotor(5,n) = m;
                            Rotor(6,n) = area_pad;
                            n = n + 1;
                        end
                    end
                end
                end
                
                u=u+1;
                theta = theta + 1;
            end
            
            Ri = Ri + .001;
        end
        c=c+1;
        thcknss = thcknss + .001;
    end
    
    
    R0 = R0 + .001;
end
    
%%
%lightest Models
   
c=1;
mass_last = 100;
for n = 1:length(Rotor(1,:))
    
    mass = (Rotor(5,n));
    
    if mass < mass_last 
        Rotor_least_weight(5,c) = Rotor(5,n);
        Rotor_least_weight(1,c) = Rotor(1,n);
        Rotor_least_weight(2,c) = Rotor(2,n);
        Rotor_least_weight(3,c) = Rotor(3,n);
        Rotor_least_weight(4,c) = Rotor(4,n);
        Rotor_least_weight(6,c) = Rotor(6,n);
        
        mass_last = mass;
    elseif mass == mass_last
        Rotor_least_weight(5,c+1) = Rotor(5,n);
        Rotor_least_weight(1,c+1) = Rotor(1,n);
        Rotor_least_weight(2,c+1) = Rotor(2,n);
        Rotor_least_weight(3,c+1) = Rotor(3,n);
        Rotor_least_weight(4,c+1) = Rotor(4,n);
        Rotor_least_weight(6,c+1) = Rotor(6,n);
        
        c=c+1;
    end
end

%%
%Fitness

fitness = 1;
highest_fitness = 0;
for n = 1:length(Rotor(1,:))
    
    R0 = Rotor(1,n);
    thcknss = Rotor(2,n);
    Ri = Rotor(3,n);
    theta = Rotor(4,n);
    mass = Rotor(5,n);
    area = Rotor(6,n);
    theta_radians = theta * pi/180;
    
    theta_fit = 100 - (abs(30 - theta) / ((30 + theta)/2)) * 100;
    Ri_fit = 100 - (abs( .18 - Ri) / (( .18 + Ri)/2) * 100);
    braking_torque_fit = .5 * area * P_max * ((R0+Ri)/2) / 100;
    area_fit = area ;
    mass_fit = 126 - (mass + (area*pad_density));
    temp_fit = -.5 * theta_radians * mew * P_max * Ri * (Ri^2 - R0^2) / 10;
    
    fitness = 1 * theta_fit + 1 * Ri_fit + 1 * area_fit + 1 * braking_torque_fit + 1 * mass_fit + 1 * temp_fit;
    
    if fitness > highest_fitness
        best(1,1) = Rotor(1,n);
        best(2,1) = Rotor(2,n);
        best(3,1) = Rotor(3,n);
        best(4,1) = Rotor(4,n);
        best(5,1) = Rotor(5,n);
        best(6,1) = Rotor(6,n);
        c=n;
        highest_fitness = fitness;
    end
end    
  
    
    
    
    
    
    
    
    
    
    
    
    