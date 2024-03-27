#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "Channel.h"
#include <boost/math/distributions/rayleigh.hpp>

//X is zero mean Gaussian distributed random variable with a standard deviation of 10.0 dB
std::normal_distribution<double> Gaussian (0.0,10);    //normal distribution 即 Gaussian distribution
boost::math::rayleigh_distribution<double> rayleigh(0.8); // standard rayleigh distribution
std::uniform_real_distribution<double> random_p(0.0, 1.0); // uniform random variable between 0.0 and 1.0 for inverse transform sampling
std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
int random_count = 100000;
double X=0.0,H=0.0,p;
bool first_time = true;

//弧度轉角度
double RadToDeg(const double & radian){
    return radian * 180 / pi;
}

//角度轉弧度
double DegToRad(const double & degree){
    return degree * pi / 180;
}

//計算整個Channel gain matrix
void Calculate_Channel_Gain_Matrix(std::vector<AP_Node> & APlist, std::vector<UE_Node> & UElist ,std::vector<std::vector<double>> & Channel_Gain_Matrix){
    for (int i = 0 ; i < AP_Num ; i++) {
        for (int j = 0 ; j < UE_Num ; j++) {
            if(i < RF_AP_Num)
                Channel_Gain_Matrix[i][j] = Estimate_one_RF_Channel_Gain(APlist[i], UElist[j]);
            else
                Channel_Gain_Matrix[i][j] = Estimate_one_VLC_Channel_Gain2(APlist[i], UElist[j]);
        }
    }
}


//計算某個 pair(RF_AP,UE)的Channel gain
/*
    RF Channel gain 公式 ：
    Gμ,α(f) = sqrt(10^(−L(d)/10)) * hr,

    hr is modelled by a Rayleigh variate on each OFDM sub-carrier with variance 2.46 dB

    L(d) = L(d0) + 10vlog10(d/d0) + X,
        L(d0) = 47.9 dB
        d0 = 1 m
        v = 1.6
        X is the shadowing component which is assumed to be a zero mean Gaussian distributed random variable with a standard deviation of 1.8 dB
 */
double Estimate_one_RF_Channel_Gain(AP_Node ap, UE_Node ue){

    double d = GetDistance_AP_UE(ap, ue);
    //calculate the X and H that are used for this simulation by averaging 100000 random results of the corresponding distribution
    if(first_time == true)
    {
        for(int i=0;i<random_count;i++)
        {
            X += Gaussian(gen);
            p = random_p(gen);
            H += quantile(rayleigh,p);
        }
        X /= random_count;
        H /= random_count;
        first_time = false;
    }
    double L_d;
    // free-space path loss function
    if(d < d_ref)
        L_d = 20 * log10(central_freq * d) - 147.5;
    else
        L_d = 20 * log10(central_freq * pow(d,2.75)/pow(d_ref,1.75)) - 147.5;

    double rf_channel_gain =  pow(H,2) * pow(10,((-1)*L_d+X)/10.0);

    return rf_channel_gain;
}

//計算某個 pair(VLC_AP,UE)的Channel gain
double Estimate_one_VLC_Channel_Gain(AP_Node ap, UE_Node ue){
    double incidence_angle = Get_Incidence_Angle_AP_UE(ap, ue);

    //若入射角 >= FOV/2則Channel gain爲0
    if(RadToDeg(incidence_angle) >= VLC_field_of_view )
        return 0.0;

    //lambertian raidation coefficient m = -ln2 / ln(cos(Φ1/2))
    const double lambertian_coefficient = (-1) / (log2(cos(DegToRad(VLC_PHI_half))));    //cos只吃弧度,所以要轉

    const double irradiance_angle = incidence_angle;

    //取得AP的位置
    Vector VLC_AP_Pos = ap.GetPosition();

    //取得UE的位置
    Vector UE_Pos = ue.GetPosition();

    const double height_diff = VLC_AP_Pos.z - UE_Pos.z;

    const double distance = GetDistance_AP_UE(ap, ue);

    double channel_gain = VLC_receiver_area * (lambertian_coefficient + 1) / (2 * pi * pow(distance,2)); // first term

    channel_gain *= pow(cos(irradiance_angle) , lambertian_coefficient); // second term

    channel_gain *= VLC_filter_gain; // third  term

    channel_gain *= pow(VLC_refractive_index , 2) / pow(sin(DegToRad(VLC_field_of_view)) , 2); // fourth term

    channel_gain *= cos(incidence_angle); // last term

    return channel_gain;
}

double Estimate_one_VLC_Channel_Gain2(AP_Node ap, UE_Node ue) {
    const double incidence_angle = Get_Incidence_Angle_AP_UE(ap, ue); // IMPORTANT: this is in RADIAN

    // no channel if exceeds FOV
    if (RadToDeg(incidence_angle) >= VLC_field_of_view) {
        return 0.0;
    }

    // TODO (alex#2#): alpha should be global constant to reduce computation.
    const double lambertian_coefficient = -(log2(2) / (log2(cos(DegToRad(VLC_PHI_half))))); //cos只吃弧度,所以要轉;

    const double irradiance_angle = incidence_angle;

    const double distance = GetDistance_AP_UE(ap, ue);
    //if (ap.GetID() == 9 && ue.GetID() == 0)
        //std::cout << distance <<std::endl;

    double channel = (lambertian_coefficient + 1) * VLC_receiver_area / (2 * pi * pow(distance, 2)); // first term

    channel = channel * pow(cos(irradiance_angle), lambertian_coefficient); // second term

    channel = channel * VLC_filter_gain * VLC_concentrator_gain; // third and fourth term

    channel = channel * cos(incidence_angle); // last term

    return channel;
}

double GetDistance_AP_UE(AP_Node ap, UE_Node ue){

    //取得AP的位置
    Vector AP_Pos = ap.GetPosition();

    //取得UE的位置
    Vector UE_Pos = ue.GetPosition();

    //高度差 h
    const double height_diff = AP_Pos.z - UE_Pos.z;

    //xy在平面上的距離 p
    double dx = AP_Pos.x - UE_Pos.x;
    double dy = AP_Pos.y - UE_Pos.y;

    const double plane_diff = sqrt(pow(dx,2)+pow(dy,2));

    //歐基米德距離
    return sqrt(pow(height_diff,2)+pow(plane_diff,2));
}

double Get_Incidence_Angle_AP_UE(Ptr<Node> AP, Ptr<Node> UE){

    //取得AP的位置
    Ptr<MobilityModel> AP_MobilityModel = AP->GetObject<MobilityModel> ();
    Vector AP_Pos = AP_MobilityModel->GetPosition ();

    //取得UE的位置
    Ptr<MobilityModel> UE_MobilityModel = UE->GetObject<MobilityModel> ();
    Vector UE_Pos = UE_MobilityModel->GetPosition ();


    //高度差 h
    const double height_diff = AP_Pos.z - UE_Pos.z;

    //xy在平面上的距離 p
    double dx = AP_Pos.x - UE_Pos.x;
    double dy = AP_Pos.y - UE_Pos.y;
    const double plane_diff = sqrt(pow(dx,2)+pow(dy,2));

    //斜邊 b
    const double hypotenuse = sqrt(pow(height_diff,2)+pow(plane_diff,2));

    //這裡採用餘弦定理來算角度，θ=cos^(-1)( (a^2+b^2-c^2)/2ab )
    const double angle = acos((pow(height_diff,2) + pow(hypotenuse,2) - pow(plane_diff,2)) / (2*height_diff*hypotenuse));


    //回傳的是角度
    return angle;
}


double Get_Incidence_Angle_AP_UE(AP_Node ap, UE_Node ue){

    //取得AP的位置
    Vector AP_Pos = ap.GetPosition();

    //取得UE的位置
    Vector UE_Pos = ue.GetPosition();


    //高度差 h
    const double height_diff = AP_Pos.z - UE_Pos.z;

    //xy在平面上的距離 p
    double dx = AP_Pos.x - UE_Pos.x;
    double dy = AP_Pos.y - UE_Pos.y;
    const double plane_diff = sqrt(pow(dx,2)+pow(dy,2));

    //斜邊 b
    const double hypotenuse = sqrt(pow(height_diff,2)+pow(plane_diff,2));

    //這裡採用餘弦定理來算角度，θ=cos^(-1)( (a^2+b^2-c^2)/2ab )
    const double angle = acos((pow(height_diff,2) + pow(hypotenuse,2) - pow(plane_diff,2)) / (2*height_diff*hypotenuse));


    //回傳的是角度
    return angle;
}
