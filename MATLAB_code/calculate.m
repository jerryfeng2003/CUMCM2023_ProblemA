function eta=calculate(eta,x,y,h_mirro,d_mirro,z0,data)
h_tower=80;
h_collector=8;
d_collector=7;
data_size=60;
h0_target=h_tower-h_collector/2;
delta_h=h_collector/4;
ang_extend=4.65e-3;
eta_ref=0.92;

cudown01=zeros(60,4);
for i = 1:data_size
    
        target=(-1)^round(rand(1,1)*10) *delta_h+h0_target;
        distance=(x^2+y^2+(target-z0)^2)^(1/2);
        eta_at=0.99321-0.0001176*distance+1.97*10^(-8)*distance^2;
        cudown01(i,1)=eta_at;
        sun_altitude=data(i,4);
        sun_azimuth=data(i,5);
        DNI=data(i,6);
        sun_vec=[sin(sun_azimuth)*cos(sun_altitude) cos(sun_azimuth)*cos(sun_altitude) sin(sun_azimuth)];
        collect_vec=-[x y z0-target]/distance;
        coll_azimuth=acos(dot([0 -1],(-[x y])/norm([x y])));
        cos_2_beta=dot(sun_vec,collect_vec);
        cos_beta=(1/2*(cos_2_beta+1))^(1/2);
        eta_cos=cos_beta;
        cudown01(i,2)=eta_cos;

        eta_sb=0.993753037;
        cudown01(i,3)=eta_sb;

        eta_truc=0.451944401;
        cudown01(i,4)=eta_truc;
          
end

result_temp=zeros(60,4);
result_temp=sum(cudown01,1)/data_size;

eta=sum(result_temp,1)/length(result_temp(:,1));