% data=xlsread('solar_data.xlsx');
% z0=4;
% z_end=6;
% data_size=60;
% test_size=20;
% h_mirro0=8;
% d_mirro0=8;
% h_mirro_start=2;
% d_mirro_start=2;
% eta_ref=0.92;
% d_pass=5;
% total=0;
% line0=10000;
% DNI_aver=sum(data(:,6),1)/data_size;
engage=1e+4;
set=400;
base=600;


location01=[0 101 z0 0];
%add function
eta=[0.9784    0.6618    1.0000    0.9557];
location01(1,4)=DNI_aver*eta_ref*prod(eta)*1e-6*h_mirro0*d_mirro0;
total=total+location01(1,4);

while(total<=61)
    fprintf('total=%d\n',total);
    test_d0=(total/61)^engage*base+set;
    i=1;
    best=zeros(test_size,4);
    x0=location01(length(location01(:,1)),1);
    y0=location01(length(location01(:,1)),2);
    x=x0;
    y=y0;
    test_set=zeros(test_size,2);
    while(i<=test_size)
        x=x0;
        y=y0;
        test_d=rand(1,1)*test_d0;
        ang=rand(1,1)*2*pi;
        while(test_d<d_mirro0+d_pass)
            test_d=rand(1,1)*test_d0;
        end
        lines=line0;
        if(length(location01(:,1))<lines)
            lines=length(location01(:,1));
        end
        while sum((sum(abs([x y].*ones(size(location01(:,1:2)))-location01(:,1:2)).^2,2).^(1/2)<=(d_mirro0+d_pass)* ...
                ones(size(location01(:,1)))))||(norm([x y])<100)||(norm([x y])>base+set)
            %||(ang<=pi/6)||(ang>=5*pi/6)
            r=test_d+d_pass;
            ang=rand(1,1)*2*pi;
            x=x0+r*cos(ang);
            y=y0+r*sin(ang);
        end
        test_set(i,:)=[x y];
        i=i+1;
    end
    for i=1:test_size
        %add function
        eta=[0 0 0 0];
        eta=calculate(eta,x,y,h_mirro0,d_mirro0,z0,data);
        best(i,1:2)=test_set(i,:);
        best(i,3)=z0;
        best(i,4)=DNI_aver*eta_ref*prod(eta)*1e-6*h_mirro0*d_mirro0;
    end
    best=sortrows(best,4,'descend');
    location01=[location01;best(1,:)];
    total=total+best(1,4);

end

location_t=location01(:,1:2);
z_result=0;
for z=z0:0.1:8
    total_temp=0;
    for i=1:length(location_t(:,1))
        %add function
        eta=[0 0 0 0];
        eta=calculate(eta,x ,y,h_mirro0,d_mirro0,z0,data);
        total_temp=total_temp+DNI_aver*eta_ref*prod(eta)*1e-6*h_mirro0*d_mirro0;
    end
    if(total_temp>total)
        total=total_temp;
        z_result=z;
    end
end

h_mirro=h_mirro0;
d_mirro=d_mirro0;
for h=h_mirro0:0.5:8
    for d=d_mirro0+1:0.5:8
        total_temp=0;
        for i=1:length(location_t(:,1))
            %add function
            eta=[0 0 0 0];
            eta=calculate(eta,x,y,h_mirro0,d_mirro0,z0,data);
            total_temp=total_temp+DNI_aver*eta_ref*prod(eta)*1e-6*h_mirro0*d_mirro0;
        end
        if(total_temp>total)
            total=total_temp;
            h_mirro=h;
            d_mirro=d;
        end
    end
end

%add function:final arguments