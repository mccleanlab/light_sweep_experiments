function Msn2 = Msn2_CT(t,Msn2_params)

signal_type = Msn2_params.signal_type;
A = Msn2_params.A; % amplitude
t0 = Msn2_params.t0; % t_start
t1 = Msn2_params.t1; % t_pulse
t2 = Msn2_params.t2; % t_interpulse
cycles = Msn2_params.cycles; % # pulses
c1 = Msn2_params.c1; % import rate
c2 = Msn2_params.c2; % export rate

Msn2 = zeros(length(t), 1);

for i=1:length(t)
    if signal_type==1 % Single pulse Msn2
        if (t(i) < t0)
            Msn2(i) = 0;
        elseif and((t(i) >= t0), (t(i) < (t0+t1)))
            Msn2(i) = A.*(1 - exp(-c1.*(t(i)-t0)));
        else
            Msn2(i) = A.*(1 - exp(-c1.*(t(i)-t0))).*exp(-c2.*(t(i) - (t0 + t1)));
        end               
    elseif signal_type==2 % Repeat pulses of Msn2
        if (t(i) < t0)
            Msn2(i) = 0;
        elseif ((t(i)-t0-floor((t(i)-t0)/(t1+t2))*(t1+t2)) <= t1) && (t(i)<(t0 + cycles*(t1+t2)))
            Msn2(i) = A.*(1 - exp(-c1.*(t(i)-t0-floor((t(i)-t0)/(t1+t2))*(t1+t2))));
        elseif ((t(i)-t0-floor((t(i)-t0)/(t1+t2))*(t1+t2)) >= t1) && (t(i)<(t0 + cycles*(t1+t2)))
            Msn2(i) = A.*(1 - exp(-c1.*(t(i)-t0-floor((t(i)-t0)/(t1+t2))*(t1+t2)))).*exp(-c2.*((t(i)-t0-floor((t(i)-t0)/(t1+t2))*(t1+t2)) - (t0 + t1)));
        else
            Msn2(i) = 0;
        end        
    elseif signal_type==3 % Step Msn2
        if (t(i) < t0)
            Msn2(i) = 0;
        elseif and((t(i) >= t0), (t(i) < (t0+t1)))
            Msn2(i) = A;
        else
            Msn2(i) = 0;
        end
    end
end

Msn2(Msn2<0) = 0;



