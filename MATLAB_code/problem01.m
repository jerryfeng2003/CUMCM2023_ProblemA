data=xlsread('solar_data.xlsx');
location01=xlsread('data.xlsx');
h_tower=80;
h_collector=8;
d_collector=7;
z0=6;
data_size=60;
h0_target=h_tower-h_collector/2;
delta_h=h_collector/4;
h_mirro=8;
d_mirro=4;
ang_extend=4.65e-3;
eta_ref=0.92;

cudown01=zeros(length(location01(:,1)),4,60);
for i = 1:data_size
    fprintf('range=%d\n',i);
    for j=1:length(location01(:,1))
%         h_mirro=location02(i,5);
%         d_mirro=location02(i,6);
        target=(-1)^j*delta_h+h0_target;
        x=location01(j,1);
        y=location01(j,2);
        locate=[x y z0];
        distance=(x^2+y^2+(target-z0)^2)^(1/2);
        eta_at=0.99321-0.0001176*distance+1.97*10^(-8)*distance^2;
        cudown01(j,1,i)=eta_at;

        sun_altitude=data(i,4);
        sun_azimuth=data(i,5);
        DNI=data(i,6);
        sun_vec=[sin(sun_azimuth)*cos(sun_altitude) cos(sun_azimuth)*cos(sun_altitude) sin(sun_azimuth)];
        collect_vec=-[x y z0-target]/distance;
        coll_azimuth=acos(dot([0 -1],(-[x y])/norm([x y])));
        cos_2_beta=dot(sun_vec,collect_vec);
        cos_beta=(1/2*(cos_2_beta+1))^(1/2);
        eta_cos=cos_beta;
        cudown01(j,2,i)=eta_cos;

        d_unsort=sum(abs(location01-[x y]).^2,2).^(1/2);
        ang_unsort=abs(acos((location01-[x y])*[0;-1]./d_unsort)-sun_azimuth);
        mat_unsort=[d_unsort ang_unsort location01];
        mat_sorted=sortrows(mat_unsort,1,"ascend");
        for k=2:length(mat_sorted(:,1))
            if(mat_sorted(k,2)<pi/4)
                cover01=[mat_sorted(k,3:4) z0];
                break;
            end
        end
        mirro0_vec=(sun_vec+collect_vec)/norm(sun_vec+collect_vec);
        cov01_coll_vec=-cover01/norm(cover01);
        mirro1_vec=(sun_vec+cov01_coll_vec)/norm(sun_vec+cov01_coll_vec);
        T1=[sin(sun_azimuth) -cos(sun_azimuth) 0;cross(sun_vec,[sin(sun_azimuth) -cos(sun_azimuth) 0]);sun_vec];
        T1=T1';
        mirro0_extend_x=cross([mirro0_vec(1:2) 0],[0 0 1])/norm([mirro0_vec(1:2) 0])*h_mirro/2;
        mirro0_A1=locate+mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A2=locate+mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A3=locate-mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A4=locate-mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro1_extend_x=cross([mirro1_vec(1:2) 0],[0 0 1])/norm([mirro1_vec(1:2) 0])*h_mirro/2;
        mirro1_A1=cover01+mirro1_extend_x+cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A2=cover01+mirro1_extend_x-cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A3=cover01-mirro1_extend_x-cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A4=cover01-mirro1_extend_x+cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro; 
        locate_sun=locate*T1;
        cover01_sun=cover01*T1;
        mirro0_A1_sun=mirro0_A1*T1;
        mirro0_A2_sun=mirro0_A2*T1;
        mirro0_A3_sun=mirro0_A3*T1;
        mirro0_A4_sun=mirro0_A4*T1;
        mirro1_A1_sun=mirro1_A1*T1;
        mirro1_A2_sun=mirro1_A2*T1;
        mirro1_A3_sun=mirro1_A3*T1;
        mirro1_A4_sun=mirro1_A4*T1;
        mirro0_A1_sun=mirro0_A1_sun(1:2);
        mirro0_A2_sun=mirro0_A2_sun(1:2);
        mirro0_A3_sun=mirro0_A3_sun(1:2);
        mirro0_A4_sun=mirro0_A4_sun(1:2);
        mirro1_A1_sun=mirro1_A1_sun(1:2);
        mirro1_A2_sun=mirro1_A2_sun(1:2);
        mirro1_A3_sun=mirro1_A3_sun(1:2);
        mirro1_A4_sun=mirro1_A4_sun(1:2);
        poly0=polyshape([mirro0_A1_sun(1) mirro0_A2_sun(1) mirro0_A3_sun(1) mirro0_A4_sun(1)],[mirro0_A1_sun(2) mirro0_A2_sun(2) mirro0_A3_sun(2) mirro0_A4_sun(2)]);
        poly1=polyshape([mirro1_A1_sun(1) mirro1_A2_sun(1) mirro1_A3_sun(1) mirro1_A4_sun(1)],[mirro1_A1_sun(2) mirro1_A2_sun(2) mirro1_A3_sun(2) mirro1_A4_sun(2)]);
        poly_out1=intersect(poly0,poly1);
        if(poly_out1.area==0)
            sb01=0;
        else
            sb01=poly_out1.area/poly0.area;
        end

        ang_unsort=abs(acos((location01-[x y])*[0;-1]./d_unsort)-coll_azimuth);
        mat_unsort=[d_unsort ang_unsort location01];
        mat_sorted=sortrows(mat_unsort,1,"ascend");
        for k=2:length(mat_sorted(:,1))
            if(mat_sorted(k,2)<pi/4)
                cover01=[mat_sorted(k,3:4) z0];
                break;
            end
        end
        cov01_coll_vec=-cover01/norm(cover01);
        mirro1_vec=(sun_vec+cov01_coll_vec)/norm(sun_vec+cov01_coll_vec);
        T1=[sin(coll_azimuth) -cos(coll_azimuth) 0;cross(collect_vec,[sin(coll_azimuth) -cos(coll_azimuth) 0]);collect_vec];
        T1=T1';
        mirro0_extend_x=cross([mirro0_vec(1:2) 0],[0 0 1])/norm([mirro0_vec(1:2) 0])*h_mirro/2;
        mirro0_A1=locate+mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A2=locate+mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A3=locate-mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A4=locate-mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro1_extend_x=cross([mirro1_vec(1:2) 0],[0 0 1])/norm([mirro1_vec(1:2) 0])*h_mirro/2;
        mirro1_A1=cover01+mirro1_extend_x+cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A2=cover01+mirro1_extend_x-cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A3=cover01-mirro1_extend_x-cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro;
        mirro1_A4=cover01-mirro1_extend_x+cross(mirro1_extend_x,mirro1_vec)/h_mirro*d_mirro; 
        locate_coll=locate*T1;
        cover01_coll=cover01*T1;
        mirro0_A1_coll=mirro0_A1*T1;
        mirro0_A2_coll=mirro0_A2*T1;
        mirro0_A3_coll=mirro0_A3*T1;
        mirro0_A4_coll=mirro0_A4*T1;
        mirro1_A1_coll=mirro1_A1*T1;
        mirro1_A2_coll=mirro1_A2*T1;
        mirro1_A3_coll=mirro1_A3*T1;
        mirro1_A4_coll=mirro1_A4*T1;
        mirro0_A1_coll=mirro0_A1_coll(1:2);
        mirro0_A2_coll=mirro0_A2_coll(1:2);
        mirro0_A3_coll=mirro0_A3_coll(1:2);
        mirro0_A4_coll=mirro0_A4_coll(1:2);
        mirro1_A1_coll=mirro1_A1_coll(1:2);
        mirro1_A2_coll=mirro1_A2_coll(1:2);
        mirro1_A3_coll=mirro1_A3_coll(1:2);
        mirro1_A4_coll=mirro1_A4_coll(1:2);
        poly0=polyshape([mirro0_A1_coll(1) mirro0_A2_coll(1) mirro0_A3_coll(1) mirro0_A4_coll(1)],[mirro0_A1_coll(2) mirro0_A2_coll(2) mirro0_A3_coll(2) mirro0_A4_coll(2)]);
        poly1=polyshape([mirro1_A1_coll(1) mirro1_A2_coll(1) mirro1_A3_coll(1) mirro1_A4_coll(1)],[mirro1_A1_coll(2) mirro1_A2_coll(2) mirro1_A3_coll(2) mirro1_A4_coll(2)]);
        poly_out2=intersect(poly0,poly1);
        if(poly_out2.area==0)
            sb02=0;
        else
            sb02=poly_out2.area/poly0.area;
        end
        
        poly0=polyshape([mirro0_A1(1) mirro0_A2(1) mirro0_A3(1) mirro0_A4(1)],[mirro0_A1(2) mirro0_A2(2) mirro0_A3(2) mirro0_A4(2)]);
        sun_vec_down=[sun_vec(1:2) 0]/norm([sun_vec(1:2) 0]);
        s1=cross(sun_vec_down,[0 0 1])*d_collector/2;
        s2=-s1;
        k=dot(sun_vec,sun_vec_down)/dot(sun_vec,[0 0 1]);
        s3=s2-sun_vec_down*h_tower*k;
        s4=-sun_vec_down*(d_collector/2+h_tower*k);
        s5=s1-sun_vec_down*h_tower*k;
        poly1=polyshape([s1(1) s2(1) s3(1) s4(1) s5(1)],[s1(2) s2(2) s3(2) s4(2) s5(2)]);
        poly_out3=intersect(poly0,poly1);
        if(poly_out3.area==0)
            sb03=0;
        else
            sb03=poly_out3.area/poly0.area;
        end       

        eta_sb=1-sb01-sb02-sb03;
        if(eta_sb<0)
            eta_sb=0;
        end
        cudown01(j,3,i)=eta_sb;

        collect_vec_down=[collect_vec(1:2) 0]/norm([collect_vec(1:2) 0]);
        mirro0_extend_x=cross([mirro0_vec(1:2) 0],[0 0 1])/norm([mirro0_vec(1:2) 0])*(h_mirro/2+distance*tan(ang_extend));
        mirro0_A1=locate+mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A2=locate+mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A3=locate-mirro0_extend_x-cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        mirro0_A4=locate-mirro0_extend_x+cross(mirro0_extend_x,mirro0_vec)/h_mirro*d_mirro;
        collect_vec_down_add=collect_vec_down*d_collector/2;
        s1=cross(collect_vec_down_add,[0 0 1])+[0 0 h_tower-h_collector];
        s2=-cross(collect_vec_down_add,[0 0 1])+[0 0 h_tower-h_collector];
        s3=-cross(collect_vec_down_add,[0 0 1])+[0 0 h_tower];
        s4=cross(collect_vec_down_add,[0 0 1])+[0 0 h_tower];   
        T2=[0 0 1;cross(collect_vec_down,[0 0 1]);collect_vec_down];
        T2=T2';
        mirro0_A1_coll=mirro0_A1*T2;
        mirro0_A2_coll=mirro0_A2*T2;
        mirro0_A3_coll=mirro0_A3*T2;
        mirro0_A4_coll=mirro0_A4*T2;
        s1_coll=s1*T2;
        s2_coll=s2*T2;
        s3_coll=s3*T2;
        s4_coll=s4*T2;
        poly0=polyshape([mirro0_A1_coll(1)+target mirro0_A2_coll(1)+target mirro0_A3_coll(1)+target mirro0_A4_coll(1)+target],[mirro0_A1_coll(2) mirro0_A2_coll(2) mirro0_A3_coll(2) mirro0_A4_coll(2)]);
        poly1=polyshape([s1_coll(1) s2_coll(1) s3_coll(1) s4_coll(1)],[s1_coll(2) s2_coll(2) s3_coll(2) s4_coll(2)]);
        poly_out4=intersect(poly0,poly1);
        if(eta_sb==0)
            eta_trunc=0;
        else
            eta_truc=poly_out4.area/(poly0.area*eta_sb);
        end
        if(eta_truc>1)
            eta_trunc=1;
        end
        cudown01(j,4,i)=eta_truc;
          
    end
end

result_temp=zeros(60,4);
for i=1:data_size
    result_temp(i,:)=sum(cudown01(:,:,i),1)/length(location01(:,1));
end

result01=zeros(12,5);
for i=0:11
    eta_temp=sum(result_temp(5*i+1:5*i+5,:),1)/5;
    DNI_temp=sum(data(5*i+1:5*i+5,6))/5;
    result01(i+1,:)=[eta_ref*prod(eta_temp) eta_temp(2) eta_temp(3) eta_temp(4) eta_ref*prod(eta_temp)*DNI_temp];
end

result02_temp=sum(result01,1)/12;
result02=[result02_temp(1:4) result02_temp(5)*h_mirro*d_mirro*length(location01(:,1)) result02_temp(5)];

result01(:,5)=result01(:,5)/1e+3;
result02(6)=result02(6)/1e+3;
result02(5)=result02(5)/1e+6;

writematrix(result01, 'problem01_table01_.xlsx');
writematrix(result02, 'problem01_table02_.xlsx');