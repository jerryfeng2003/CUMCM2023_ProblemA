G0 = 1366; 
latitude = 39.4; 
longitude = 98.5; 
altitude = 3000; 
H = altitude / 1000; 

time_points = [9, 10.5, 12, 13.5, 15];

data = zeros(60, 6); 

for month = 1:12
    for i = 1:length(time_points)
        D = (month - 1) * 30 + 21;
        
        ST = time_points(i);
        omega = pi / 12 * (ST - 12);
        
        delta = asin(sin(2 * pi * D / 365) * sin(2 * pi / 365 * 23.45));
        
        alpha_s = asin(cos(delta) * cosd(latitude) * cos(omega) + sind(delta) * sind(latitude));
        gamma_s = acos((sin(delta) - sin(alpha_s) * sind(latitude)) / (cos(alpha_s) * cosd(latitude)));
        
        a = 0.4237 - 0.00821 * (6 - H)^2;
        b = 0.5055 + 0.00595 * (6.5 - H)^2;
        c = 0.2711 + 0.01858 * (2.5 - H)^2;
        DNI = G0 * (a + b * exp(-c / sin(alpha_s)));
        
        row = (month - 1) * 5 + i;
        data(row, 1) = month;
        data(row, 2) = D;
        data(row, 3) = ST;
        data(row, 4) = alpha_s;
        data(row, 5) = gamma_s;
        data(row, 6) = DNI;
     
    end
end

data_table = array2table(data, 'VariableNames', {'Month', 'Day', 'SolarTime', 'SolarAltitude', 'SolarAzimuth', 'DNI'});

writetable(data_table, 'solar_data.xlsx');
