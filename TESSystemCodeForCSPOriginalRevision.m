clc
clear all
format short g
%% Start Code
%Input Parameters
% Thermophysical Properties of Paraffin Wax as PCM with Geometric of TES
pi=3.14;%Pi Number
H=0.608;%1.1;% Height of tank (m)
D=0.380;%0.9;% Diameter oftank (m)
e=0.47;%0.379;% Porosity
d_p=0.076;%0.05;% Diameter of PCM capsule
mdat=0.0796;%0.0416; Mass Flow Rate for Fluid (Changeable);(NOTICE:0.0416 for 150 & 0.0796 for 200)
ro_s=670;%838;% PCM density in solid mode
ro_l=640;%834;% PCM density in liquid mode
C_ps=2400;% PCM specific heat in solid mode
C_pl=1800;% PCM specific heat in liquid mode
k_s=0.4;% PCM thermal conductivity in solid mode
k_l=0.2;% PCM thermal conductivity in liquid mode
mu_pcm=5e-6;% PCM dynamic viscosity
h_fgpcm=142.7;%254;% PCM Latent heat
h_fgpcmSolid=38;
T_m1=319.15;%326.05;% Peak temperature of the PCM during the solid-solid transition
T_m2=333.15;% Peak temperature of the PCM during the solid-solid transition
fi=input('Volume fraction of nano particle:');%Volume Fraction of Nano Particle
if fi>0
    program1=1;%input('picking out of type of nanoparticles (1=Au; 2=Ag; 3=Ni; 4=TiO2; 5=Al):');
    p1=program1;
    if p1==1
        ro_p=19300;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=310;%Thermal Conductivity of Nano Particle
        C_p_p=129;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==2
        ro_p=10.50e3;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=419;%Thermal Conductivity of Nano Particle
        C_p_p=235;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==3
        ro_p=8900;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=77.5e-9;%Diameter of Nano Particle
        k_p=90.7;%Thermal Conductivity of Nano Particle
        C_p_p=444;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==4
        ro_p=4230;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=8.5;%Thermal Conductivity of Nano Particle
        C_p_p=690;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==5
        ro_p=2702;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=237;%Thermal Conductivity of Nano Particle
        C_p_p=903;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    end
end
%% Start of calculations
T_0Ch=389;%input('Value of input temperature for heat transfer fluid:');% input temperature for heat transfer fluid (K)
T_pcm_in=318.27;%input('Value of input temperature for PCM:');% input temperature for PCM (K)
e2=H/d_p;% Number of cells
e3=(D/2)/(d_p/2);% Number of radius elements
T_ifN11=T_0Ch;
T_isr1=T_pcm_in;
T_ifN1_outi1=0;
T_ifN_outi1N1=0;
T_ifN_outi11=zeros(e2,e3);% (output)
T_isr_outi11=zeros(e2,e3);% (output)
for j=1:e2% N elements
    R=D/2;
    r_0=R/e3;
    r_p=r_0;
    for k=1:e3% Radial elements
        if T_isr1<T_m1
            ro_pcm1=ro_s;
            C_ppcm1=C_ps;
            k_pcm1=k_s;
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr1;
                T_isr_outi1r1=0;
            end
            mu_f1=((2.1897e-11)*(T_ifN11.^4))-((3.055e-8)*(T_ifN11.^3))+((1.6028e-5)*(T_ifN11.^2))-(0.0037524*T_ifN11)+0.33158;%Fluid Dynamic Viscosity
            C_p_f1=((1.1105e-5)*(T_ifN11.^3))-(0.0031078*(T_ifN11.^2))-(1.478*T_ifN11)+4631.9;%Fluid Specific Heat
            k_f1=((1.5362e-8)*(T_ifN11.^3))-((2.261e-05)*(T_ifN11.^2))+(0.010879*T_ifN11)-1.0294;%Fluid Thermal Conductivity
            ro_f1=((-1.5629e-5)*(T_ifN11.^3))+(0.011778*(T_ifN11.^2))-(3.0726*T_ifN11)+1227.8;%Fluid Density
            A_i=pi*(r_p.^2);
            u_f1=mdat/(A_i*ro_f1);
            Pr_f1=(C_p_f1*mu_f1)/k_f1;
            Re_f1=(ro_f1*u_f1*(d_p))/(mu_f1);
            h_f1=(2+(1.1*((6*(1-e)).^0.6)*(Re_f1.^0.6)*(Pr_f1.^(1/3))))*k_pcm1/d_p;%(2+(0.6*(Re_f1.^(1/2))*(Pr_f1.^(1/3))))*(k_pcm1/d_p);%
            A11=(ro_f1*C_p_f1*e);
            A21=(ro_f1*C_p_f1*u_f1)/(d_p);
            A31=(k_f1*e)/(d_p.^2);
            A41=((6*(1-e)*h_f1)/d_p);
            A51=ro_pcm1*C_ppcm1;
            A61=(k_pcm1/(r_p.^2));
            A71=(6*h_f1)/d_p;
            B11=A11*T_ifN11;
            B21=(-A21+A31)*T_ifN1_outi1;
            B31=(A21+A31)*T_ifN_outi1N1;
            B41=A41*T_isr1;
            B51=A61*T_isr_outi1r1;
            B61=A61*T_isr1_outi1;
            B71=A51*T_isr1;
            B81=A71*T_ifN11;
            A111=A11+(2*A31)+A41;
            B111=B11+B21+B31+B41;
            A221=A51+(2*A61)+A71;
            B221=B51+B61+B71+B81;
            %X=T_ifN_outi1,T_isr_outi3
            A1111=[A111,0;0,A221];
            B1111=[B111;B221];
            X111=A1111\B1111;
            T_ifN_outi11(j,k)=X111(1,:);
            T_isr_outi11(j,k)=X111(2,:);
            T_isr1=T_isr_outi11(j,k);
            T_ifN11=T_ifN_outi11(j,k);
            Capacity1(j,k)=C_p_f1;
            Cp1=Capacity1;
            HTC1(j,k)=h_f1;
            HTC1Ch=HTC1;
            if T_isr1>T_m1
                T_isr2=T_isr1;
                T_ifN21=T_ifN11;
            end
        elseif T_isr2>T_m1 && T_isr2<T_m2
            mu_f2=((2.1897e-11)*(T_ifN21.^4))-((3.055e-8)*(T_ifN21.^3))+((1.6028e-5)*(T_ifN21.^2))-(0.0037524*T_ifN21)+0.33158;%Fluid Dynamic Viscosity
            C_p_f2=((1.1105e-5)*(T_ifN21.^3))-(0.0031078*(T_ifN21.^2))-(1.478*T_ifN21)+4631.9;%Fluid Specific Heat
            k_f2=((1.5362e-8)*(T_ifN21.^3))-((2.261e-05)*(T_ifN21.^2))+(0.010879*T_ifN21)-1.0294;%Fluid Thermal Conductivity
            ro_f2=((-1.5629e-5)*(T_ifN21.^3))+(0.011778*(T_ifN21.^2))-(3.0726*T_ifN21)+1227.8;%Fluid Density
            if fi==0
                k_pcm2=(k_s+k_l)/2;
            else
                k_nf2=k_f2*((k_p+((n-1)*k_f2)-((n-1)*fi*(k_f2-k_p)))/(k_p+((n-1)*k_f2)+(fi*(k_f2-k_p))));
                k_pcm2=(k_s+k_l+k_nf2)/2;
            end
            A_i=pi*(r_p.^2);
            u_f2=mdat/(A_i*ro_f2);
            Pr_f2=(C_p_f2*mu_f2)/k_f2;
            Re_f2=(ro_f2*u_f2*(d_p))/(mu_f2);
            h_f2=(2+(1.1*((6*(1-e)).^0.6)*(Re_f2.^0.6)*(Pr_f2.^(1/3))))*k_pcm2/d_p;%(2+(0.6*(Re_f2.^(1/2))*(Pr_f2.^(1/3))))*(k_pcm2/d_p);%
            A12=(ro_f2*C_p_f2*e);
            A22=(ro_f2*C_p_f2*u_f2)/(d_p);
            A32=(k_f2*e)/(d_p.^2);
            if fi==0
                ro_pcm2=ro_s;
                h_fg=h_fgpcm;
                C_ppcm2=((C_ps+C_pl)/2)+(h_fg/(T_m2-T_m1));
            else
                ro_nf2=(ro_p*fi)+(ro_f2*(1-fi));
                ro_pcm2=ro_s+ro_nf2;
                h_fg=(ro_f2*h_fgpcm*(1-fi))/ro_nf2;
                C_ppcm2=((C_ps+C_pl)/2)+(h_fg/(T_m2-T_m1));
            end
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr2;
                T_isr_outi1r1=0;
            end
            A42=((6*(1-e)*h_f2)/d_p);
            A52=ro_pcm2*C_ppcm2;
            A62=(k_pcm2/(r_p.^2));
            A72=(6*h_f2)/d_p;
            B12=A12*T_ifN21;
            B22=(-A22+A32)*T_ifN1_outi1;
            B32=(A22+A32)*T_ifN_outi1N1;
            B42=A42*T_isr2;
            B52=A62*T_isr_outi1r1;
            B62=A62*T_isr1_outi1;
            B72=A52*T_isr2;
            B82=A72*T_ifN21;
            A112=A12+(2*A32)+A42;
            B112=B12+B22+B32+B42;
            A222=A52+(2*A62)+A72;
            B222=B52+B62+B72+B82;
            %X=T_ifN_outi1,T_isr_outi1
            A1112=[A112,0;0,A222];
            B1112=[B112;B222];
            X112=A1112\B1112;
            T_ifN_outi11(j,k)=X112(1,:);
            T_isr_outi11(j,k)=X112(2,:);
            T_isr2=T_isr_outi11(j,k);
            T_ifN21=T_ifN_outi11(j,k);
            Capacity2(j,k)=C_p_f2;
            Cp2=Capacity2;
            HTC2(j,k)=h_f2;
            HTC2Ch=HTC2;
            if T_isr2>T_m2
                T_isr3=T_isr2;
                T_ifN31=T_ifN21;
            end
        elseif T_isr3>T_m2
            ro_pcm3=ro_l;
            C_ppcm3=C_pl;
            k_pcm3=k_l;
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr3;
                T_isr_outi1r1=0;
            end
            mu_f3=((2.1897e-11)*(T_ifN31.^4))-((3.055e-8)*(T_ifN31.^3))+((1.6028e-5)*(T_ifN31.^2))-(0.0037524*T_ifN31)+0.33158;%Fluid Dynamic Viscosity
            C_p_f3=((1.1105e-5)*(T_ifN31.^3))-(0.0031078*(T_ifN31.^2))-(1.478*T_ifN31)+4631.9;%Fluid Specific Heat
            k_f3=((1.5362e-8)*(T_ifN31.^3))-((2.261e-05)*(T_ifN31.^2))+(0.010879*T_ifN31)-1.0294;%Fluid Thermal Conductivity
            ro_f3=((-1.5629e-5)*(T_ifN31.^3))+(0.011778*(T_ifN31.^2))-(3.0726*T_ifN31)+1227.8;%Fluid Density
            A_i=pi*(r_p.^2);
            u_f3=mdat/(A_i*ro_f3);
            Pr_f3=(C_p_f3*mu_f3)/k_f3;
            Re_f3=(ro_f3*u_f3*(d_p))/(mu_f3);
            h_f3=(2+(1.1*((6*(1-e)).^0.6)*(Re_f3.^0.6)*(Pr_f3.^(1/3))))*k_pcm3/d_p;%(2+(0.6*(Re_f3.^(1/2))*(Pr_f3.^(1/3))))*(k_pcm3/d_p);%
            A13=(ro_f3*C_p_f3*e);
            A23=(ro_f3*C_p_f3*u_f3)/(d_p);
            A33=(k_f3*e)/(d_p.^2);
            A43=((6*(1-e)*h_f3)/d_p);
            A53=ro_pcm3*C_ppcm3;
            A63=(k_pcm3/(r_p.^2));
            A73=(6*h_f3)/d_p;
            B13=A13*T_ifN31;
            B23=(-A23+A33)*T_ifN1_outi1;
            B33=(A23+A33)*T_ifN_outi1N1;
            B43=A43*T_isr3;
            B53=A63*T_isr_outi1r1;
            B63=A63*T_isr1_outi1;
            B73=A53*T_isr3;
            B83=A73*T_ifN31;
            A113=A13+(2*A33)+A43;
            B113=B13+B23+B33+B43;
            A223=A53+(2*A63)+A73;
            B223=B53+B63+B73+B83;
            %X=T_ifN_outi1,T_isr_outi1
            A1113=[A113,0;0,A223];
            B1113=[B113;B223];
            X113=A1113\B1113;
            T_ifN_outi11(j,k)=X113(1,:);
            T_isr_outi11(j,k)=X113(2,:);
            T_isr3=T_isr_outi11(j,k);
            T_ifN31=T_ifN_outi11(j,k);
            Capacity3(j,k)=C_p_f3;
            Cp3=Capacity3;
            HTC3(j,k)=h_f3;
            HTC3Ch=HTC3;
        end
        r_p=r_p+r_0;
    end
end
%CpChTotal=[Cp1;Cp2;Cp3];
%CpChargingSUM=sum(CpChTotal(:));
T_f_out_Ch=T_ifN_outi11(e2,e3);
T_s_out_Ch=T_isr_outi11(e2,e3);
Q_Ch_HTF=mdat*((Capacity1(1,1)*T_0Ch)-(Capacity3(e2,e3)*T_f_out_Ch));%mdat*CpChargingSUM*(T_0Ch-T_f_out_Ch);%
Q_Ch_PCM=(mdat*C_ppcm1*(T_isr2-T_pcm_in))+(mdat*C_ppcm2*(T_isr3-T_isr2))+(mdat*C_ppcm3*(T_s_out_Ch-T_isr3))+(mdat*h_fgpcm)+(mdat+h_fgpcmSolid)
disp 'The energy stored in storage in Charging mode'
Q_Storage_Ch=Q_Ch_HTF+Q_Ch_PCM
%% Calculations pertinent to Discharging
% Thermophysical Properties of Paraffin Wax as PCM with Geometric of TES
pi=3.14;%Pi Number
H=0.608;%1.1;% Height of tank (m)
D=0.380;%0.9;% Diameter oftank (m)
e=0.47;%0.379;% Porosity
d_p=0.076;%0.05;% Diameter of PCM capsule
mdat=0.0796;%0.0416; Mass Flow Rate for Fluid (Changeable);(NOTICE:0.0416 for 150 & 0.0796 for 200)
ro_s=670;%838;% PCM density in solid mode
ro_l=640;%834;% PCM density in liquid mode
C_ps=2400;% PCM specific heat in solid mode
C_pl=1800;% PCM specific heat in liquid mode
k_s=0.4;% PCM thermal conductivity in solid mode
k_l=0.2;% PCM thermal conductivity in liquid mode
mu_pcm=5e-6;% PCM dynamic viscosity
h_fgpcm=142.7;%254;% PCM Latent heat
h_fgpcmSolid=38;
T_m1=319.15;%326.05;% Peak temperature of the PCM during the solid-solid transition
T_m2=333.15;% Peak temperature of the PCM during the solid-solid transition
fi=input('Volume fraction of nano particle:');%Volume Fraction of Nano Particle
if fi>0
    program1=1;%input('picking out of type of nanoparticles (1=Au; 2=Ag; 3=Ni; 4=TiO2; 5=Al):');
    p1=program1;
    if p1==1
        ro_p=19300;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=310;%Thermal Conductivity of Nano Particle
        C_p_p=129;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==2
        ro_p=10.50e3;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=419;%Thermal Conductivity of Nano Particle
        C_p_p=235;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==3
        ro_p=8900;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=77.5e-9;%Diameter of Nano Particle
        k_p=90.7;%Thermal Conductivity of Nano Particle
        C_p_p=444;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==4
        ro_p=4230;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=8.5;%Thermal Conductivity of Nano Particle
        C_p_p=690;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    elseif p1==5
        ro_p=2702;%Density of Nano Particle
        ro_f_ref=999;%Fluid Density in the 20 degree celcius
        d_pn=67.5e-9;%Diameter of Nano Particle
        k_p=237;%Thermal Conductivity of Nano Particle
        C_p_p=903;%Specific Heat for Nano Paricle
        n=3;%empirical shape factor of Nano Particle
        Av=6.022*(10e23);%Avogadro
        M=0.0180153;%Mole Mass of base Fluid
    end
end
%% Start of calculations
T_0Dis=T_f_out_Ch;%input('Value of input temperature for heat transfer fluid:');% input temperature for heat transfer fluid (K)
T_pcm_in=T_s_out_Ch;%input('Value of input temperature for PCM:');% input temperature for PCM (K)
e2=H/d_p;% Number of cells
e3=(D/2)/(d_p/2);% Number of radius elements
T_ifN12=T_0Dis;
T_isr1=T_pcm_in;
T_ifN1_outi1=0;
T_ifN_outi1N1=0;
T_ifN_outi12=zeros(e2,e3);% (output)
T_isr_outi12=zeros(e2,e3);% (output)
for j=1:e2% N elements
    R=D/2;
    r_0=R/e3;
    r_p=r_0;
    for k=1:e3% Radial elements
        if T_isr1>T_m2
            ro_pcm1=ro_l;
            C_ppcm1=C_pl;
            k_pcm1=k_l;
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr1;
                T_isr_outi1r1=0;
            end
            mu_f1=((2.1897e-11)*(T_ifN12.^4))-((3.055e-8)*(T_ifN12.^3))+((1.6028e-5)*(T_ifN12.^2))-(0.0037524*T_ifN12)+0.33158;%Fluid Dynamic Viscosity
            C_p_f1=((1.1105e-5)*(T_ifN12.^3))-(0.0031078*(T_ifN12.^2))-(1.478*T_ifN12)+4631.9;%Fluid Specific Heat
            k_f1=((1.5362e-8)*(T_ifN12.^3))-((2.261e-05)*(T_ifN12.^2))+(0.010879*T_ifN12)-1.0294;%Fluid Thermal Conductivity
            ro_f1=((-1.5629e-5)*(T_ifN12.^3))+(0.011778*(T_ifN12.^2))-(3.0726*T_ifN12)+1227.8;%Fluid Density
            A_i=pi*(r_p.^2);
            u_f1=mdat/(A_i*ro_f1);
            Pr_f1=(C_p_f1*mu_f1)/k_f1;
            Re_f1=(ro_f1*u_f1*(d_p))/(mu_f1);
            h_f1=(2+(1.1*((6*(1-e)).^0.6)*(Re_f1.^0.6)*(Pr_f1.^(1/3))))*k_pcm1/d_p;%(2+(0.6*(Re_f1.^(1/2))*(Pr_f1.^(1/3))))*(k_pcm1/d_p);%
            A11=(ro_f1*C_p_f1*e);
            A21=(ro_f1*C_p_f1*u_f1)/(d_p);
            A31=(k_f1*e)/(d_p.^2);
            A41=((6*(1-e)*h_f1)/d_p);
            A51=ro_pcm1*C_ppcm1;
            A61=(k_pcm1/(r_p.^2));
            A71=(6*h_f1)/d_p;
            B11=A11*T_ifN12;
            B21=(-A21+A31)*T_ifN1_outi1;
            B31=(A21+A31)*T_ifN_outi1N1;
            B41=-A41*T_isr1;
            B51=A61*T_isr_outi1r1;
            B61=A61*T_isr1_outi1;
            B71=A51*T_isr1;
            B81=-A71*T_ifN12;
            A111=A11+(2*A31)-A41;
            B111=B11+B21+B31+B41;
            A221=A51+(2*A61)-A71;
            B221=B51+B61+B71+B81;
            %X=T_ifN_outi1,T_isr_outi1
            A1111=[A111,0;0,A221];
            B1111=[B111;B221];
            X111=A1111\B1111;
            T_ifN_outi12(j,k)=X111(1,:);
            T_isr_outi12(j,k)=X111(2,:);
            T_isr1=T_isr_outi12(j,k);
            T_ifN12=T_ifN_outi12(j,k);
            Capacity1(j,k)=C_p_f1;
            Cp1=Capacity1;
            if T_isr1<T_m2
                T_isr2=T_isr1;
                T_ifN22=T_ifN12;
            end
        elseif T_isr2>T_m1 && T_isr2<T_m2
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr2;
                T_isr_outi1r1=0;
            end
            mu_f2=((2.1897e-11)*(T_ifN22.^4))-((3.055e-8)*(T_ifN22.^3))+((1.6028e-5)*(T_ifN22.^2))-(0.0037524*T_ifN22)+0.33158;%Fluid Dynamic Viscosity
            C_p_f2=((1.1105e-5)*(T_ifN22.^3))-(0.0031078*(T_ifN22.^2))-(1.478*T_ifN22)+4631.9;%Fluid Specific Heat
            k_f2=((1.5362e-8)*(T_ifN22.^3))-((2.261e-05)*(T_ifN22.^2))+(0.010879*T_ifN22)-1.0294;%Fluid Thermal Conductivity
            ro_f2=((-1.5629e-5)*(T_ifN22.^3))+(0.011778*(T_ifN22.^2))-(3.0726*T_ifN22)+1227.8;%Fluid Density
            if fi==0
                k_pcm2=(k_s+k_l)/2;
            else
                k_nf2=k_f2*((k_p+((n-1)*k_f2)-((n-1)*fi*(k_f2-k_p)))/(k_p+((n-1)*k_f2)+(fi*(k_f2-k_p))));
                k_pcm2=(k_s+k_l+k_nf2)/2;
            end
            A_i=pi*(r_p.^2);
            u_f2=mdat/(A_i*ro_f2);
            Pr_f2=(C_p_f2*mu_f2)/k_f2;
            Re_f2=(ro_f2*u_f2*(d_p))/(mu_f2);
            h_f2=(2+(1.1*((6*(1-e)).^0.6)*(Re_f2.^0.6)*(Pr_f2.^(1/3))))*k_pcm2/d_p;%(2+(0.6*(Re_f2.^(1/2))*(Pr_f2.^(1/3))))*(k_pcm2/d_p);%
            A12=(ro_f2*C_p_f2*e);
            A22=(ro_f2*C_p_f2*u_f2)/(d_p);
            A32=(k_f2*e)/(d_p.^2);
            if fi==0
                ro_pcm2=ro_s;
                h_fg=h_fgpcm;
                C_ppcm2=((C_ps+C_pl)/2)+(h_fg/(T_m2-T_m1));
            else
                ro_nf2=(ro_p*fi)+(ro_f2*(1-fi));
                ro_pcm2=ro_s+ro_nf2;
                h_fg=(ro_f2*h_fgpcm*(1-fi))/ro_nf2;
                C_ppcm2=((C_ps+C_pl)/2)+(h_fg/(T_m2-T_m1));
            end
            A42=((6*(1-e)*h_f2)/d_p);
            A52=ro_pcm2*C_ppcm2;
            A62=(k_pcm2/(r_p.^2));
            A72=(6*h_f2)/d_p;
            B12=A12*T_ifN22;
            B22=(-A22+A32)*T_ifN1_outi1;
            B32=(A22+A32)*T_ifN_outi1N1;
            B42=-A42*T_isr2;
            B52=A62*T_isr_outi1r1;
            B62=A62*T_isr1_outi1;
            B72=A52*T_isr2;
            B82=-A72*T_ifN22;
            A112=A12+(2*A32)-A42;
            B112=B12+B22+B32+B42;
            A222=A52+(2*A62)-A72;
            B222=B52+B62+B72+B82;
            %X=T_ifN_outi1,T_isr_outi1
            A1112=[A112,0;0,A222];
            B1112=[B112;B222];
            X112=A1112\B1112;
            T_ifN_outi12(j,k)=X112(1,:);
            T_isr_outi12(j,k)=X112(2,:);
            T_isr2=T_isr_outi12(j,k);
            T_ifN22=T_ifN_outi12(j,k);
            Capacity2(j,k)=C_p_f2;
            Cp2=Capacity2;
            if T_isr2<T_m1
                T_isr3=T_isr2;
                T_ifN32=T_ifN22;
            end
        elseif T_isr3<T_m1
            ro_pcm3=ro_s;
            C_ppcm3=C_ps;
            k_pcm3=k_s;
            if k==1
                T_isr1_outi1=0;
                T_isr_outi1r1=0;
            else
                T_isr1_outi1=T_isr3;
                T_isr_outi1r1=0;
            end
            mu_f3=((2.1897e-11)*(T_ifN32.^4))-((3.055e-8)*(T_ifN32.^3))+((1.6028e-5)*(T_ifN32.^2))-(0.0037524*T_ifN32)+0.33158;%Fluid Dynamic Viscosity
            C_p_f3=((1.1105e-5)*(T_ifN32.^3))-(0.0031078*(T_ifN32.^2))-(1.478*T_ifN32)+4631.9;%Fluid Specific Heat
            k_f3=((1.5362e-8)*(T_ifN32.^3))-((2.261e-05)*(T_ifN32.^2))+(0.010879*T_ifN32)-1.0294;%Fluid Thermal Conductivity
            ro_f3=((-1.5629e-5)*(T_ifN32.^3))+(0.011778*(T_ifN32.^2))-(3.0726*T_ifN32)+1227.8;%Fluid Density
            A_i=pi*(r_p.^2);
            u_f3=mdat/(A_i*ro_f3);
            Pr_f3=(C_p_f3*mu_f3)/k_f3;
            Re_f3=(ro_f3*u_f3*(d_p))/(mu_f3);
            h_f3=(2+(1.1*((6*(1-e)).^0.6)*(Re_f3.^0.6)*(Pr_f3.^(1/3))))*k_pcm3/d_p;%(2+(0.6*(Re_f3.^(1/2))*(Pr_f3.^(1/3))))*(k_pcm3/d_p);%
            A13=(ro_f3*C_p_f3*e);
            A23=(ro_f3*C_p_f3*u_f3)/(d_p);
            A33=(k_f3*e)/(d_p.^2);
            A43=((6*(1-e)*h_f3)/d_p);
            A53=ro_pcm3*C_ppcm3;
            A63=(k_pcm3/(r_p.^2));
            A73=(6*h_f3)/d_p;
            B13=A13*T_ifN32;
            B23=(-A23+A33)*T_ifN1_outi1;
            B33=(A23+A33)*T_ifN_outi1N1;
            B43=-A43*T_isr3;
            B53=A63*T_isr_outi1r1;
            B63=A63*T_isr1_outi1;
            B73=A53*T_isr3;
            B83=-A73*T_ifN32;
            A113=A13+(2*A33)-A43;
            B113=B13+B23+B33+B43;
            A223=A53+(2*A63)-A73;
            B223=B53+B63+B73+B83;
            %X=T_ifN_outi1,T_isr_outi1
            A1113=[A113,0;0,A223];
            B1113=[B113;B223];
            X113=A1113\B1113;
            T_ifN_outi12(j,k)=X113(1,:);
            T_isr_outi12(j,k)=X113(2,:);
            T_isr3=T_isr_outi12(j,k);
            T_ifN32=T_ifN_outi12(j,k);
            Capacity3(j,k)=C_p_f3;
            Cp3=Capacity3;
        end
        r_p=r_p+r_0;
    end
end
T_f_out_Disch=T_ifN_outi12(e2,e3);
T_s_out_Disch=T_isr_outi12(e2,e3);
%CpDischTotal=[Cp1;Cp2;Cp3];
%CpDischargingSUM=sum(CpDischTotal(:));
Q_Disch_HTF=mdat*((Capacity3(e2,e3)*T_f_out_Disch)-(Capacity3(e2,e3)*T_0Dis));%mdat*CpDischargingSUM*(T_f_out_Disch-T_0Dis);%
Q_Disch_PCM=(mdat*C_ppcm1*(T_pcm_in-T_isr2))+(mdat*C_ppcm2*(T_isr2-T_isr3))+(mdat*C_ppcm3*(T_isr3-T_s_out_Disch))+(mdat*h_fgpcm)+(mdat*h_fgpcmSolid);
disp 'The efficiency of storage:'
E=Q_Disch_HTF/Q_Ch_HTF