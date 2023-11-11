data=xlsread('solar_data.xlsx');
z0=4.2;
z_end=6;
data_size=60;
test_size=20;
h_mirro0=8;
d_mirro0=8;
h_mirro_start=2;
d_mirro_start=2;
eta_ref=0.92;
d_pass=5;
total=0;
line0=10000;
DNI_aver=sum(data(:,6),1)/data_size;
engage=1e+1;
set=300;
base=700;


location01=[0 101 z0 0 h_mirro0 d_mirro0];
%add function
eta=[0.9784    0.6618    1.0000    0.9557];
location01(1,4)=DNI_aver*eta_ref*prod(eta)*1e-6*h_mirro0*d_mirro0;
total=total+location01(1,4);

while(total<=62)
    fprintf('total=%d\n',total);
    test_d0=(total/62)^engage*base+set;
    i=1;
    best=ones(test_size,6);
    best(:,5)=best(:,5)*h_mirro0;
    best(:,6)=best(:,6)*d_mirro0;
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
    temp_result=best(1,4);
    temp=0;
    for z=z0:0.1:z_end
        for h=h_mirro_start:0.1:h_mirro0
            for d=d_mirro_start:0.5:d_mirro0
                eta=calculate(eta,best(1,1),best(1,2),h,d,z,data);
                temp=DNI_aver*eta_ref*prod(eta)*h*1e-6*d;
                if(temp/(h*d)>temp_result/(best(1,5)*best(1,6)))
                    tmep_result=temp;
                    best(1,3)=z;
                    best(1,5)=h;
                    best(1,6)=d;
                end
            end
        end
    end
    best(1,4)=temp_result;
    location01=[location01;best(1,:)];  
    total=total+best(1,4);

end
scatter(location01(:,1),location01(:,2));