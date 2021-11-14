function dict=epg_fisp_mrf(RFpulses ,TR_times ,TE, T1, T2)
%simulation of FISP using EPG. All times are in seconds and angles are in
%Radians

Ntr = length(TR_times);	% Number of TRs
Nstates = 20;	% Number of states to simulate 

P = zeros(3,Nstates);	% State matrix
P(3,1)=-1;		% Equilibrium magnetization.

dict=zeros(1,length(TR_times)); % Vector holding the recieved signals

for k = 1:Ntr
    TR=TR_times(k);
    flipang=abs(RFpulses(k));
    flipphase=angle(RFpulses(k));
    
    P = epg_rf(P,flipang,flipphase);		% RF pulse
    
    % FID for time TE and then sample signal
    P = epg_grelax(P,T1,T2,TE,0,0,0,0);   % Relaxation.  这次脉冲激发后产生的信号
    dict(k) = P(1,1);      			% Signal is F0 state.
     
    % -- Simulate relaxation and spoiler gradient
    P = epg_grelax(P,T1,T2,TR-TE,1,0,1,1);   % spoiler gradient, relaxation.  %下一次脉冲前，信号的初始状态
    
end;



