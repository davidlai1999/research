#ifndef GLOBAL_ENVIROMENT_H
#define GLOBAL_ENVIROMENT_H

#define DEBUG_MODE 0    //1=show debug message, 0=don't show
#define DEBUG_DETAIL 0  //more debug messages,can also show the movement of the UE , only support with 1 UE. (haven't use this for a while, don't know if this still works)
#define BENCHMARK  0    //0=proposed, 1=benchmark
#define NO_DEMAND 0 //0=normal environment, 1=environment without users' required rates

//the minimum rate that is set for the inactive links, unit: bps
const uint64_t min_rate = 100000;

const double pi = 3.14159265359;

////////////////////////////////////////////////////////
/////////        Simulation Constants          /////////
////////////////////////////////////////////////////////
//simulation time, unit: seconds
const int simulation_time = 300.0; // 300

//time between each execution of the algorithm, unit: seconds
const double delta_t = 0.05;

//average vertical handover delay, unit: seconds
const double VHO_delay = 0.5;

//average horizontal handover delay, unit: seconds
const double HHO_delay = 0.2;

//time to trigger, unit: seconds
const double time_TTT = 0.16;

//handover margin, unit: db
const int delta_HOM = 1;

//the parameter of the benchmark that is adjusted according to the avg_speed, if avg_speed<=1.0 then lambda = 1.0
//, if 1.0<avg_speed<1.6 then lambda = avg_speed*2-1.0, if avg_speed=1.6 then lambda = 2.4, if avg_speed=1.7 then lambda = 4.0, else lambda = 10.0
const double lambda = 1.0;

//average speed of  UEs, unit: m/s
const double avg_speed =1.0;

//the maximum pause time for the RWP model, unit: seconds
const double pause_time = 0.0;

//LS(lowest satisfaction) of the proposed algorithm
const double min_guaranteed_portion = 0.7;

//the step size of the adjustment of the dynamically adjusted satisfaction level of UEs
const double delta_p = 0.05;

//this proportion of the UEs are bounded to an area
const double crowded_portion = 0;

////////////////////////////////////////////////////////
/////////        ENVIRONMENT CONSTANTS          ////////
////////////////////////////////////////////////////////

// room size(for both width and length)
const double room_size=10;

// - RF AP number
const int RF_AP_Num = 4;

// - VLC AP number
const int VLC_AP_Num = 16;

// - Total AP number
const int AP_Num = RF_AP_Num + VLC_AP_Num;

// - UE number
const int UE_Num = 20;

// VLC AP y axis start
const double VLC_AP_Y = 3.75;

// VLC AP x axis start
const double VLC_AP_X = 3.75;

// VLC AP SPACE
const double VLC_AP_SPACE = 1.25;

// RF AP y axis start
const double RF_AP_Y = 2.5;

// RF AP x axis start
const double RF_AP_X = 2.5;

// RF AP SPACE
const double RF_AP_SPACE = 5.0;



////////////////////////////////////////////////////////
/////////          Each RF AP                  /////////
////////////////////////////////////////////////////////

// - height
const double RF_AP_height = 3;

// - power (17 dBm ~= 0.05W)
const double RF_AP_Power = 0.05;

// - RF bandwidth
const int RF_AP_Bandwidth = 10;

// -174dBm/Hz
// - power spectral density
const double RF_AWGN_spectral_density = 1e-16;

// reference distance for free-space path loss function
const double d_ref = 5.0;

// central carrier frequency
const double central_freq = 5.0e9;

// - Max power
const double RF_Max_Power = 10.0;

////////////////////////////////////////////////////////
/////////          Each VLC AP                  ////////
////////////////////////////////////////////////////////
// - Max power
const double VLC_Max_Power = 10.0;

// - height
const double VLC_AP_height = 3;

// - modulated optical power
const double VLC_AP_Popt = 1;

// - VLC bandwidth
const int VLC_AP_Bandwidth = 75;

// - AWGN power spectral density    10^-22  A^2/Hz = 10^-16  A^2/MHz
const double  VLC_AWGN_spectral_density = 1e-16;

// - Detector responsivity = 0.53 A/W
const double kappa = 0.53 ;

const int reuse_factor = 2; // RB id range from 0 to g_frequency_reuse_factor-1

const double bandwidth_per_cell = VLC_AP_Bandwidth / reuse_factor;

////////////////////////////////////////////////////////
/////////          Each UE                      ////////
////////////////////////////////////////////////////////

//   - height
const double UE_height = 0.0;

//   minimum required data rate
const double min_required_rate = 1;//Mbps

//   maximum required data rate
const double max_required_rate = 150;//Mbps

////////////////////////////////////////////////////////
/////////         VLC CHANNEL CONSTANTS         ////////
////////////////////////////////////////////////////////

// semi-angle of FOV
const double VLC_field_of_view = 50;

// semi-angle at half-illumination (phi_1/2) 60
const int VLC_PHI_half = 60;

// gain of optical filter (g_of(psi))         1
const int VLC_filter_gain = 1;

// gain of optical concentrator (g_oc(psi))   1, can use this for simplicity, or use refractive index to calculate the actual value
const int VLC_concentrator_gain = 1;
// refractive index
const double VLC_refractive_index = 1.5;

// physical area for PD receiver              1 cm^2 = 0.0001 m^2
const double VLC_receiver_area = 0.0001;

// reflection efficiency (rho)                0.75
const double VLC_reflect_efficiency = 0.75;

const int number_of_layer = 4;

// BER_Lower_Bound
const double BER_LB = 0.00001;

const double speed_threshold = 1.4;

#endif
