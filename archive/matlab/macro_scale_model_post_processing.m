%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% macro_scale_model_post_processing.m
%% Dr. Brittany Bannish, University of Central Oklahoma
%% Stephen B. Gregg, University of Central Oklahoma
%% Updated on 2/8/2017
%% 
%% This MATLAB file calculates and plots results from a macro scale model
%% of blood clot degredation. 
%% 
%% Output:
%%  
%% i. Calculations
%%      (1) front_velocity
%%      (2) std_front_velocity
%%      (3) degradation_rate
%%      (4) std_degradation_rate
%%      
%% ii. Plots
%%      (1) Fiber Height Degradation: degradation height vs time
%%      (2) Fiber Degradation: fraction of fibers degraded vs time
%%      (3) Total Fibers Degraded: number of fibers degraded vs time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
format long

% binary formatting: type sun for sun and lin for linux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  arch = 'lin';
  if arch == 'sun',
     binaryformat = 'ieee-be';
  else 
     binaryformat = 'ieee-le';    
  end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% global variables
%%
%% N: number of nodes in one row in horizontal direction 
%% F: number of nodes in one column in vertical direction
%% number_of_simulations: number of independent macroscale model simulati-
%%   ons performed
%% tf: final time in macroscale model (in seconds)
%% degradable_fibers_constant: number of degradable fibers in the clot
%%
%% Note: the above quantities (N, F, number_of_simulations, final_simulation_time) should match the correspo-
%% nding numbers from the macroscale model code that was used to generate
%% the data
%%
%% total_fibers: total number of fibers in clot (including ghost edges)
%%               (this counts the number of fibers in the lattice diagram)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=342;%93;%87;%93;%117;%93;%83;
F=39;%121;%113;%121;%152;%121;%108;
number_of_simulations=10;
final_simulation_time=15*60;%20*60;
number_of_degradable_fibers_constant=28700;%25761;%40833;%25761;%20501; %this is num-enoFB from macro Fortran code

total_fibers=(2*N-1)*F+N*(F-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% open files created with Fortran
%%
%% fids(k) '[file]'_[tPA conc]_[PLG conc]_[insert]_[insert]_[quartile]
%%
%% fids(1) 'lastmove_': gives the last time step at which the front moves for each x-location
%% fids(2) 'Nsave_': number of time points per simulation
%% fids(7) 'deg_': 
%% fids(8) 'move': 
%% fids(9) 'plot':
%% fids(10) 'tsave': 
%%
%% number_of_time_points: gives the total number of entries in time-dependent vectors. So if final_simulation_time=20*60 and we save every 10 seconds, then there are 1200/10+1=121 saved time points
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fids(1) = fopen('lastmove_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
  last_move = fread(fids(1),[N,number_of_simulations],'int');
  fclose(fids(1));
fids(2) = fopen('Nsave_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
  Nsave = fread(fids(2),inf,'int');
  number_of_time_points=Nsave(1)+1;
  fclose(fids(2));
fids(7) = fopen('deg_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
fids(8) = fopen('move_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
fids(9) = fopen('plot_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
fids(10) = fopen('tsave_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);

% calculate front velocity and mean degradation rate
for ii=1:number_of_simulations
 ii
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% calculate front velocity
 %%
 %% move: gives the saved-time-step at which the front moves for each x location
 %% degradation_height: gives in each row the successive y-positions (in microns) of x-location corresponding to row number
 %% stopwatch_time: gives the simulation time at which each entry was saved
 %% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 move = fread(fids(8),[N,N],'int');
 degradation_height=fread(fids(9),[N,N],'double');
 stopwatch_time(ii,1:number_of_time_points)=fread(fids(10),number_of_time_points,'double');
  
 xy_slope=zeros(1,N);

 %to get the time in minutes, divide stopwatch_time by 60
 simulation_stopwatch_time=stopwatch_time(ii,1:number_of_time_points)/60;

 %fit a line to each x-location data (micron vs. minutes)
 for i=1:N
    %form a vector of timesteps
    timestep_of_moves=0;
    y_change_stopwatch_time=0;
    %indexmove=find(move(:,i)==0,1,'first'); %if for some reason last_move isn't correct, this seems to work....
    %timestep_of_moves=move(1:indexmove-1,i); %only look at the entries of "move" that actually correspond to the front moving. i.e., ignore 0 entries
    timestep_of_moves=move(1:last_move(i,ii),i); %only look at the entries of "move" that actually correspond to the front moving. i.e., ignore 0 entries
    for j=1:length(timestep_of_moves)
        y_change_stopwatch_time(j)=simulation_stopwatch_time(timestep_of_moves(j)); %take the times (in minutes) corresponding to the times that the front moved
    end
    
%     plot(y_change_stopwatch_time,degradation_height(i,1:indexmove-1)) %plot y-position (in microns) of the front in the ith x-location vs. the time (in minutes)
%     [linefit_slope_yint_velocity]=polyfit(y_change_stopwatch_time(2:end),degradation_height(i,2:indexmove-1),1); %fit a line to the above plot excluding the first time point
    plot(y_change_stopwatch_time,degradation_height(i,1:last_move(i,ii))) %plot y-position (in microns) of the front in the ith x-location vs. the time (in minutes)
    [linefit_slope_yint_velocity]=polyfit(y_change_stopwatch_time(2:end),degradation_height(i,2:last_move(i,ii)),1); %fit a line to the above plot excluding the first time point
    xy_slope(i)=linefit_slope_yint_velocity(1); %the slope of the line excluding the first time point gives the front velocity in microns/min

    hold all

 end
 xlabel('time (min)','FontSize',20,'FontWeight','b')
 ylabel('y-position (microns)','FontSize',20,'FontWeight','b')
 title('Fiber Height Degradation', 'FontSize', 14)
 set(gca,'FontSize',18,'FontWeight','bold')
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% front velocity calculations (excludes first time point)
 %%
 %% simulation_front_velocity: average front velocity (in microns/min) over all x-locations
 %% simulation_std_front_velocity: standard deviation across all x-locations
 %% aggregate_simulation_front_velocity(ii): stores simulation_front_velocity for each simulation
 %% aggregate_simulation_std_front_velocity(ii): stores simulation_std_front_velocity for each simulation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 simulation_front_velocity=mean(xy_slope);
 simulation_std_front_velocity=std(xy_slope);
 aggregate_simulation_front_velocity(ii)=simulation_front_velocity;
 aggregate_simulation_std_front_velocity(ii)=simulation_std_front_velocity;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS IS THE MEAN DEGRADATION RATE PART
%%
%% degradation_status: 
%%      -1: ghost edge
%%      <-1: time (in seconds) at which the fiber degraded 
%%      0: fiber is undegraded
%% 
%%
%% 
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 degradation_status_vector=fread(fids(7),number_of_time_points*total_fibers,'double');
 degradation_status_matrix=zeros(total_fibers,number_of_time_points+1);
   
 for i=2:(number_of_time_points+1)         
     degradation_status_matrix(:,i) = degradation_status_vector(((i-2)*total_fibers+1):(i-1)*total_fibers)'; %makes a big matrix with degradation states. the first column is the initial condition, and then each subsequent column is the degradation status of each edge at each saved time point
 end
     
 first_degradable_fiber=find(degradation_status_matrix(:,2)==0,1,'first'); %find the first place in the column that the degradation state is 0. i.e., the first place degradation hasn't happened. This will help identify the total number of degradable fibers in the clot
 number_of_degradable_fibers_calculated=total_fibers-first_degradable_fiber+1; %total number of degradable edges = total number of edges in clot (including ghost edges) minus the number from above plus 1. This should equal fibnum. 
 for j=2:(number_of_time_points+1)
     number_of_degraded_fibers(j-1)=length(find(degradation_status_matrix(first_degradable_fiber:end,j)<0)); %calculate the total number of degraded edges at each saved time step
 end
    
 fraction_of_degraded_fibers=number_of_degraded_fibers/number_of_degradable_fibers_calculated; %calculate the fraction of degraded edges at each saved time point
          
 aggregate_fraction_of_degraded_fibers(:,ii)=fraction_of_degraded_fibers; %matrix of fraction of degraded edges at each time
 aggregate_number_of_degraded_fibers(:,ii)=number_of_degraded_fibers; %matrix of number of degraded edges at each time
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%
%% front_velocity: average over all simulations and all x locations for each trial (micron/min; excludes first time points)
%% 
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


front_velocity=mean(aggregate_simulation_front_velocity)
std_front_velocity=mean(aggregate_simulation_std_front_velocity) 
save -ascii mspeed_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat front_velocity

%DEGRADATION RATE
figure
for j=1:number_of_simulations
    plot(stopwatch_time(j,:),aggregate_fraction_of_degraded_fibers(:,j)) %plot the fraction of fibers degraded at each time point
    hold all
end

xlabel('time (sec)','FontSize',20,'FontWeight','b')
ylabel('fraction of fibers degraded','FontSize',20,'FontWeight','b')
title('Fiber Degradation', 'FontSize', 14)
set(gca,'FontSize',18,'FontWeight','bold')

fclose(fids(7));

%find the time point at which lysis begins and ends for each run
for j=1:number_of_simulations
    timepoint_end_of_degradation=find(aggregate_fraction_of_degraded_fibers(:,j)==1,1,'first'); %Find the first place where degradation is 100% complete
    timepoint_start_of_degradation(j)=find(aggregate_fraction_of_degraded_fibers(:,j)>0,1,'first'); %Find the first place degradation starts
    if length(timepoint_end_of_degradation)==0 %if 100% lysis didn't happen in the given run
        timepoint_furthest_degradation(j)=length(aggregate_fraction_of_degraded_fibers(:,j)); %take the "end of degradation" to be the last entry of the Per_mat matrix
    else %if 100% lysis DID happen in the given run
        timepoint_furthest_degradation(j)=timepoint_end_of_degradation; %take the "end of degradation" to be the entry at which 100% lysis was first achieved
    end
    

 %define the "time" that we're interested in to be the entries/saved time
 %points for which degradation actually occurred. So start the degradation
 %time at the 1st place degradation occurred, and finish when you reach
 %100% lysis, or the end of the run, whichever comes first.
 timepoints_of_degradation=timepoint_start_of_degradation(j):timepoint_furthest_degradation(j); 

%fit a line from the first place the slope between 2 consecutive seconds is
%>=tol to the last place the slope is >=tol
mslope(1)=aggregate_fraction_of_degraded_fibers(timepoints_of_degradation(1),j); %calculate the slope between the first place degradation occurs and the time point right before that
for i=2:length(timepoints_of_degradation)
    mslope(i)=aggregate_fraction_of_degraded_fibers(timepoints_of_degradation(i),j)-aggregate_fraction_of_degraded_fibers(timepoints_of_degradation(i-1),j); %calculate the slope between all pairs of successive time points
end
%set the tolerance that we will use to determine the degradation rate.
tol=1e-3; %potentially need to make this 1e-4. fiddle with it if results look bad
left=find(mslope>=tol,1,'first'); %find the first place the slope is >= to the tolerance
right=find(mslope>=tol,1,'last'); %find the last place the slope is >= to the tolerance
if right>length(timepoints_of_degradation)
   right=length(timepoints_of_degradation);
end
timepoints_of_degradation_shifted=timepoints_of_degradation(left:right); %redefine the time so that we're only looking at the places where the slope is greater than the tolerance
[linefit_slope_yint_degradation1]=polyfit(stopwatch_time(j,timepoints_of_degradation_shifted)',aggregate_fraction_of_degraded_fibers(timepoints_of_degradation_shifted,j),1); %fit a line to the data fraction of fibers degraded vs. time (in sec)
degradation1_slope(j)=linefit_slope_yint_degradation1(1); %save this slope for each of the 10 runs
plot(stopwatch_time(j,timepoints_of_degradation_shifted)',linefit_slope_yint_degradation1(1)*(stopwatch_time(j,timepoints_of_degradation_shifted))+linefit_slope_yint_degradation1(2),'b') %plot the fitted lines

[linefit_slope_yint_degradation2]=polyfit(stopwatch_time(j,timepoints_of_degradation_shifted)'/60,aggregate_fraction_of_degraded_fibers(timepoints_of_degradation_shifted,j),1); %fit a line to the data fraction of fibers degraded vs. time (in min)
degradation2_slope(j)=linefit_slope_yint_degradation2(1); %save this slope for each of the 10 runs
end


degradation_rate=100*mean(degradation2_slope) %mean degradation rate over 10 runs, measured in percent of total fibers/min 
std_degradation_rate=100*std(degradation2_slope) %standard deviation over 10 runs

save -ascii mdeg_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat degradation_rate

%plot the number of fibers in the clot as a function of time
figure
for j=1:number_of_simulations
    plot(stopwatch_time(j,:)',number_of_degradable_fibers_constant-aggregate_number_of_degraded_fibers(:,j))
    hold all
    xlabel('time (sec)','FontSize',20,'FontWeight','b')
    ylabel('number of fibers','FontSize',20,'FontWeight','b')
    title('Total Fibers Degraded', 'FontSize', 14)
    set(gca,'FontSize',18,'FontWeight','bold')
end
save -ascii Numfib_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat aggregate_number_of_degraded_fibers

%%%%%%%For Colin:
fraction_of_undegraded_fibers=1-aggregate_fraction_of_degraded_fibers; %fraction of fibers still to be degraded in the clot
time_points=29.4*stopwatch_time(1,:)'; %the time, in seconds, scaled by 29.4 to match the size of Colin's clot

save -ascii time_points_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat time_points
save -ascii frac_undeg_fib_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat fraction_of_undegraded_fibers

fclose(fids(8));
fclose(fids(9));
fclose(fids(10));

%%%%%%%%%%Calculate mean first passage time%%%%%%%%%%%%%%%%%%%
fids(5) = fopen('mfpt_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat','r',binaryformat);
  firstpasstime = fread(fids(5),inf,'double');
fclose(fids(5)); 

%only use tPA molecules that actually hit the back of the clot
ind_pass=find(firstpasstime>0);
num_molecules_diffused_to_end=length(ind_pass) %this is the number of molecules that made it to the end
frac_molecules_diffused_to_end=length(ind_pass)/length(firstpasstime)

nonzerofpt=0;
for i=1:length(ind_pass)
    nonzerofpt(i)=firstpasstime(ind_pass(i));
end

%calculat mean first passage time by taking the mean of all the first
%passage times for the molecules that made it to the end
mean_fpt_sec=mean(nonzerofpt);
mean_fpt_min=mean_fpt_sec/60
std_min=std(nonzerofpt/60)

save -ascii meanfirst_PLG2_tPA01_finenoRBC_lessFB_F39_Q2.dat mean_fpt_min
