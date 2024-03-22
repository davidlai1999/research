#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono> // seed
#include <unistd.h>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "Algo.h"
#include "math.c" // erfinv

int cnt = 1;
int total = 0;

std::vector<double> BER_required_SINR(number_of_layer+1,0.0);

//channel gain
std::vector<std::vector<double>> Channel_Gain_Matrix(AP_Num,std::vector<double> (UE_Num,0));

//frequency reuse
std::vector<int> frequency (VLC_AP_Num,0);

//dwell time timer
std::vector<double> Dwell_timer (UE_Num, 0.0);

//timer of handover, first is VHO, second is HH0
std::vector<std::pair<double, double>> Handover_timer (UE_Num, {0.0, 0.0});
std::vector<std::pair<double, double>> Record_handover_time (UE_Num, {0.0, 0.0});

// mode 0 : UE is normal, mode 1 : UE redundant, mode 2 : UE dwell time, mode 3 : during handover
std::vector<int> UE_timer_mode (UE_Num,0);
//modulation mode
std::vector<std::vector<int>> mode_table = {
    {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 0, 3}, {0, 0, 1, 0}, {0, 0, 2, 0},
    {0, 0, 3, 0}, {0, 0, 1, 1}, {0, 1, 0, 0}, {0, 2, 0, 0}, {0, 1, 0, 1},
    {0, 1, 0, 2}, {0, 1, 0, 3}, {0, 1, 1, 0}, {1, 0, 0, 0}, {1, 0, 0, 1},
    {2, 0, 0, 0}, {1, 0, 0, 2}, {1, 0, 0, 3}, {1, 0, 1, 0}, {1, 0, 2, 0},
    {1, 0, 3, 0}, {1, 1, 0, 0}, {1, 2, 0, 0}, {2, 1, 0, 0}, {1, 1, 0, 2},
    {1, 1, 0, 3}, {2, 2, 0, 0}, {1, 1, 1, 0}, {1, 1, 2, 0}, {1, 2, 1, 0},
    {1, 2, 2, 0}, {2, 1, 1, 0}, {2, 1, 2, 0}, {2, 2, 1, 0}, {2, 2, 2, 0},
    {3, 1, 1, 0}, {2, 1, 1, 1}, {2, 2, 3, 0}, {2, 1, 1, 2}, {3, 1, 2, 0},
    {2, 1, 2, 1}, {3, 2, 1, 0}, {2, 2, 1, 1}, {3, 2, 2, 0}, {2, 2, 2, 1},
    {3, 1, 1, 1}, {3, 2, 3, 0}, {2, 2, 2, 2}, {2, 3, 1, 1}, {3, 1, 1, 2},
    {3, 1, 2, 1}, {4, 2, 1, 0}, {2, 3, 2, 1}, {3, 2, 1, 1}, {4, 2, 2, 0},
    {3, 2, 1, 2}, {3, 2, 2, 1}, {3, 3, 1, 1}, {3, 2, 3, 1}, {4, 1, 2, 1},
    {3, 3, 1, 2}, {3, 3, 2, 1}, {4, 2, 1, 1}, {3, 4, 1, 1}, {3, 3, 2, 2},
    {4, 2, 1, 2}, {4, 2, 2, 1}, {3, 4, 2, 1}, {4, 3, 1, 1}, {4, 2, 3, 1},
    {4, 3, 1, 2}, {4, 3, 2, 1}, {4, 4, 1, 1}, {4, 3, 2, 2}, {4, 3, 3, 1},
    {4, 4, 1, 2}, {4, 4, 2, 1}, {4, 3, 4, 1}, {4, 4, 2, 2}, {4, 4, 3, 1},
    {4, 4, 2, 3}, {4, 4, 3, 2}, {4, 4, 4, 1}, {4, 4, 3, 3}, {4, 4, 4, 2},
    {4, 4, 3, 4}, {4, 4, 4, 3}
};

std::vector<double> mode_rate_table;

//lists to store the APs and the UEs
std::vector<AP_Node> APlist;
std::vector<UE_Node> UElist;
std::vector<double> UE_average_rate (UE_Num,0.0);// average data rate of each UE

// data used for analysis
double simulation_steps = simulation_time / delta_t;
int UE_connect_to_AP = 0;        // The number of UE connect to AP
int Redundant_UE_count = 0;      // The number of UE redundant
int UE_connect_to_RF_AP = 0;     // The number of UE connect to RF AP
int UE_connect_to_VLC_AP = 0;    // The number of UE connect to VLC AP
int RFtoVLC = 0;                 // handover count from RF to VLC
int VLCtoRF = 0;                 // handover count from VLC to RF
int VLCtoVLC = 0;                // handover count from VLC to VLC
int RFtoRF = 0;                  // handover count from RF to RF
int VHO_count = 0;               // vertical handover count
int HHO_count = 0;               // horizontal handover count
double total_handover_delay = 0; // total time that is consumed by the handover
double total_HHO_delay = 0;      // total time that is consumed by the HHO
double total_VHO_delay = 0;      // total time that is consumed by the VHO

bool compare_channel_gain(const UE_Node& node1, const UE_Node& node2)
{
    return (Channel_Gain_Matrix[node1.Get_Associated_AP()][node1.GetID()] > Channel_Gain_Matrix[node2.Get_Associated_AP()][node2.GetID()]);
}

bool compare_residual_power(std::pair<int,double> p1, std::pair<int,double> p2)
{
    return(p1.second > p2.second);
}

// 初始化 UE 串列，並賦予每個 UE required rate
void Initialize_UE_Node_list(NodeContainer &UE_Nodes){
    for (int i = 0 ; i < UE_Num ; i++) {
        // 利用 C++ 自帶的 uniform distribution generator 生成隨機 required data rate
        std::uniform_real_distribution<double> uniform_dis(min_required_rate , max_required_rate);
        std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
        double required_data_rate = uniform_dis(re);

        //新增My_UE_Node加入myUElist中
        UElist.push_back(UE_Node(i, required_data_rate, UE_Nodes.Get(i)));
    }
}

// 初始化 AP 串列，並將相鄰 UE 的 bandwidth 錯開
void Initialize_AP_Node_list(NodeContainer & AP_Nodes){
    // Add AP into the list
    for(int i = 0 ; i < AP_Num ; i++){
        Ptr<MobilityModel> AP_MobilityModel = (AP_Nodes.Get(i))->GetObject<MobilityModel>();
        Vector pos = AP_MobilityModel->GetPosition ();
        if (i < RF_AP_Num)
            APlist.push_back(AP_Node(i, RF_Max_Power, pos, AP_Nodes.Get(i)));
        else
            APlist.push_back(AP_Node(i, VLC_Max_Power, pos, AP_Nodes.Get(i)));
    }

    //frequency reuse setting
    frequency[1]=frequency[7]=frequency[9]=frequency[15]=1;
    frequency[2]=frequency[4]=frequency[10]=frequency[12]=2;
    frequency[3]=frequency[5]=frequency[11]=frequency[13]=3;
}

void Update_UE_Condiction() {
    for (int i = 0 ; i < UElist.size() ; i++)
        UElist[i].UpdateCondiction();
}

double set_BER_Required_SINR(const int& number_of_layer) {
    /*for (int mod = 1 ; mod <= number_of_layer ; mod++) {
        double SINR = BER_LB;
        SINR = SINR * mod * pow(2, mod);
        SINR = SINR / (pow(2, mod) - 1);
        SINR = 1 - SINR;
        SINR = erfinv(SINR);
        SINR = pow(SINR, 2);
        SINR = SINR / 3 * 2 * (pow(2, 2 * mod) - 1);

        BER_required_SINR[mod-1] = (SINR);
    }*/
    BER_required_SINR[0] = 0.0;
    BER_required_SINR[1] = 9.6;
    BER_required_SINR[2] = 13.4;
    BER_required_SINR[3] = 17.8;
    BER_required_SINR[4] = 22.5;
}

double get_layer_rate(int layer, double SINR) {
    return bandwidth_per_cell / pow(2, layer) * log2(1 + SINR);
}

double get_achievable_data_rate(std::vector<int> mode) {
    double ans = 0.0;
    for (int i = 0 ; i < mode.size() ; i++) {
        if (mode[i] == 0) continue;
        ans += get_layer_rate(i+1, BER_required_SINR[mode[i]]);
    }

    return ans;
}

void create_modulation_mode_rate_table() {
    for (int i = 0 ; i < mode_table.size() ; i++) {
        double rate = get_achievable_data_rate(mode_table[i]);
        mode_rate_table.push_back(rate);
    }
}

int select_modulation_mode(double rate, int l, int r) {
    if (r - l  == 1)
        return r;

    if (l > r)
        return -1;

    int mid = (l + r)/2 ;
    if (mode_rate_table[mid] > rate)
        return select_modulation_mode(rate, l, mid);
    else if (rate > mode_rate_table[mid])
        return select_modulation_mode(rate, mid, r);
    else
        return mid;
}

void Set_UE_default_modulation_mode() {
    // select modulation mode
    for (int ue_id = 0 ; ue_id < UElist.size() ; ue_id++) {
        int mode = select_modulation_mode(UElist[ue_id].Get_Required_DataRate(), 0, mode_rate_table.size()-1);
        UElist[ue_id].Set_Mode_Index(mode);
        UElist[ue_id].Set_Modulation_Mod(mode_table[mode]);
    }
}
/**
* UE needs THIS much power from AP to use m in layer l alone
* \return (double) required power use m in l
*/
double get_layer_required_power(const int ap_id, const int ue_id,
                          const int mod, const int layer,
                          const double sum_prev_ue_power){
    /* SINR of this layer = single layer SINR using BPSK*/
    double SINR = BER_required_SINR[mod]; // REMEMBER ARRAY STARTS AT 0!!
    /* beta_l (subcarrier ratio) of this layer */
    double beta_l = (double) 1 / std::pow(2, layer);
    /* interference from prev UEs (1 to k-1) */
    double intra_I = sum_prev_ue_power * std::pow(Channel_Gain_Matrix[ap_id][ue_id], 2);
    /* AWGN */
    // double AWGN = g_total_bandwidth * g_N_0;
    double AWGN = bandwidth_per_cell * VLC_AWGN_spectral_density;
    /* inter-cell-interference I_n,k */
    double ICI = 0.0;
    /*for(int ap_b_id = 0; ap_b_id < VLC_AP_Num; ap_b_id++) {
        if (ap_b_id == ap_id) continue;
        else {
            check if same RB
            int ap_rb_id = Room::transmitter_[ap_id]->get_rb_id();
            int ap_b_rb_id = Room::transmitter_[ap_b_id]->get_rb_id();
            if (ap_rb_id == ap_b_rb_id) /* ICI ONLY if use same RB
                ICI += Room::transmitter_[ap_b_id]->get_ICI(ue_id);
        }
    }*/

    /* power for this layer */
    double power = SINR * beta_l * (intra_I + AWGN + ICI) / std::pow(Channel_Gain_Matrix[ap_id][ue_id],2);

    return power;
}

double get_VLC_minimum_required_power(const int ap_id, const int ue_id, const double sum_prev_ue_power){
    double layer_power = 0.0, first_term = 0.0, second_term = 0.0;
    for (int layer = 1 ; layer <= number_of_layer ; layer++) {
        layer_power = get_layer_required_power(ap_id, ue_id, UElist[ue_id].Get_Modulation_Mod()[layer-1], layer, sum_prev_ue_power);
        first_term += layer_power;
        second_term += sqrt(layer_power);
    }

    first_term = first_term * (pi - 1) / pi;

    second_term = pow(second_term, 2) / pi;

    double power = first_term + second_term;
    return power;
}

double get_VLC_minimum_required_power(const int ap_id, const int ue_id, std::vector<int> modulation_mode, const double sum_prev_ue_power){
    double layer_power = 0.0, first_term = 0.0, second_term = 0.0;
    for (int layer = 1 ; layer <= number_of_layer ; layer++) {
        layer_power = get_layer_required_power(ap_id, ue_id, modulation_mode[layer-1], layer, sum_prev_ue_power);
        first_term += layer_power;
        second_term += sqrt(layer_power);
    }

    first_term = first_term * (pi - 1) / pi;

    second_term = pow(second_term, 2) / pi;

    double power = first_term + second_term;
    return power;
}

int select_modulation_mode_base_on_residual_power(int ap_id, int ue_id, double residual_power, int l, int r, double sum_prev_ue_power) {
    if (r - l  == 1) {
        double l_power = get_VLC_minimum_required_power(ap_id, ue_id, mode_table[l], sum_prev_ue_power);
        return (l_power < residual_power) ? l : -1;
    }

    if (l > r)
        return -1;

    int mid = (l + r) / 2 ;
    if (get_VLC_minimum_required_power(ap_id, ue_id, mode_table[mid], sum_prev_ue_power) > residual_power) {
        return select_modulation_mode_base_on_residual_power(ap_id, ue_id, residual_power, l, mid, sum_prev_ue_power);
    }
    else if (residual_power > get_VLC_minimum_required_power(ap_id, ue_id, mode_table[mid], sum_prev_ue_power)){
        return select_modulation_mode_base_on_residual_power(ap_id, ue_id, residual_power, mid, r, sum_prev_ue_power);
    }
    else {
        return mid;
    }

}

/*int select_modulation_mode_base_on_residual_power(int ap_id, int ue_id, double residual_power, int l, int r, double sum_prev_ue_power) {
    if (r - l  == 1) {
        double l_power = get_VLC_minimum_required_power(ap_id, ue_id, mode_table[l], sum_prev_ue_power);
        double r_power = get_VLC_minimum_required_power(ap_id, ue_id, mode_table[r], sum_prev_ue_power);
        return (r_power < residual_power) ? r : ( l_power < residual_power) ? l : -1;
    }

    if (r - l  == 0)
        return l;

    if (l > r)
        return -1;

    int mid = (l + r) / 2 ;
    if (get_VLC_minimum_required_power(ap_id, ue_id, mode_table[mid], sum_prev_ue_power) > residual_power) {
        return select_modulation_mode_base_on_residual_power(ap_id, ue_id, residual_power, l, mid-1, sum_prev_ue_power);
    }
    else if (residual_power > get_VLC_minimum_required_power(ap_id, ue_id, mode_table[mid], sum_prev_ue_power)){
        return select_modulation_mode_base_on_residual_power(ap_id, ue_id, residual_power, mid, r, sum_prev_ue_power);
    }
    else {
        return mid;
    }

}*/

double get_RF_minimum_required_power(const int ap_id, const int ue_id, const double sum_prev_ue_power) {
    double channel_gain = std::pow(Channel_Gain_Matrix[ap_id][ue_id], 2);

    double result = channel_gain * sum_prev_ue_power;

    result = result + RF_AP_Bandwidth * RF_AWGN_spectral_density;

    result = result * (pow(2, UElist[ue_id].Get_Required_DataRate() / RF_AP_Bandwidth) - 1);

    result = result / channel_gain;
}

double get_RF_SINR_based_on_Residual_Power(int ap_id, int ue_id, double sum_prev_ue_power, double residual_power) {
    double channel_gain = pow(Channel_Gain_Matrix[ap_id][ue_id], 2);

    double SINR = channel_gain * residual_power;

    double AWGN = RF_AP_Bandwidth * RF_AWGN_spectral_density;

    double ICI = 0.0;

    SINR = SINR / (channel_gain * sum_prev_ue_power + AWGN + ICI);
}

double get_RF_Achievable_Data_Rate(double SINR) {
    return RF_AP_Bandwidth * log2(1 + SINR) ;
}

std::vector<double> power_allocation_after_handover(int ap_id, int ue_id) {
    if (APlist[ap_id].Get_Residual_Power() <= 0) return {0.0, APlist[ap_id].Get_Residual_Power()};
    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();

    double sum = 0.0, prev = 0.0, required_power = 0.0;
    int bigger_channel_ue = -1;
    for (int i = 0 ; i < UElist_of_AP.size() ; i++) {
        if (Channel_Gain_Matrix[ap_id][UElist_of_AP[i]] >= Channel_Gain_Matrix[ap_id][ue_id]) {
            prev = APlist[ap_id].Get_UE_Power(UElist_of_AP[i]);
            sum = sum + prev;
        }
        else {
            bigger_channel_ue = i;
            break;
        }
    }

    if (ap_id < RF_AP_Num)
        required_power = get_RF_minimum_required_power(ap_id, ue_id, sum);
    else
        required_power = get_VLC_minimum_required_power(ap_id, ue_id, sum);

    required_power = (prev > required_power) ? prev : required_power;
    double target_required_power = required_power;
    sum = sum + required_power;
    prev = required_power;

    if (bigger_channel_ue != -1) {
        for (int i = bigger_channel_ue ; i < UElist_of_AP.size() ; i++) {
            if (ap_id < RF_AP_Num) {
                required_power = get_RF_minimum_required_power(ap_id, UElist_of_AP[i], sum);
                required_power = (prev > required_power) ? prev : required_power;
                sum = sum + required_power;
                prev = required_power;
            }
            else {
                required_power = get_VLC_minimum_required_power(ap_id, UElist_of_AP[i], sum);
                required_power = (prev > required_power) ? prev : required_power;
                sum = sum + required_power;
                prev = required_power;
            }

        }
    }

    double residual_power = VLC_Max_Power - sum;
    return {target_required_power, residual_power};
}

int find_Best_VLC_Channel_Gain(int ue_id) {
    int result = -1;
    double maxi = 1e-9;
    for (int ap_id = RF_AP_Num ; ap_id < AP_Num ; ap_id++) {
        if (maxi < Channel_Gain_Matrix[ap_id][ue_id]) {
            maxi = Channel_Gain_Matrix[ap_id][ue_id];
            result = ap_id;
        }
    }

    return result;
}

int find_Best_RF_Channel_Gain(int ue_id) {
    int result = -1;
    double maxi = 1e-9;
    for (int ap_id = 0 ; ap_id < RF_AP_Num ; ap_id++) {
        if (maxi < Channel_Gain_Matrix[ap_id][ue_id]) {
            maxi = Channel_Gain_Matrix[ap_id][ue_id];
            result = ap_id;
        }
    }

    return result;
}

void Set_UE_Associated_AP_Based_On_Channel() {
    for (int i = 0 ; i < UElist.size() ; i++) {
        int ap_id = find_Best_VLC_Channel_Gain(i);

        if (ap_id == -1)
            ap_id = find_Best_RF_Channel_Gain(i);

        UElist[i].Set_Associated_AP(ap_id);
    }
}

void Update_VLC_AP_power_allocation(int ap_id) {
    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();
    std::vector<UE_Node> tmp_UElist;
    for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
        tmp_UElist.push_back(UElist[UElist_of_AP[j]]);
    }

    APlist[ap_id].Clear_Associated_UE();
    sort(tmp_UElist.begin(), tmp_UElist.end(), compare_channel_gain);

    double sum = 0.0, prev = 0.0, required_power = 0.0;
    int curr = 0;
    for (; curr < tmp_UElist.size() ; curr++) {
        int ue_id = tmp_UElist[curr].GetID();

        if (Channel_Gain_Matrix[ap_id][ue_id] == 0)
            break;

        required_power = get_VLC_minimum_required_power(ap_id, ue_id, sum);
        required_power = (prev > required_power) ? prev : required_power;
        sum = sum + required_power;

        if (sum > VLC_Max_Power) {
            sum = sum - required_power;
            if (VLC_Max_Power - sum < prev) { // condiction 1: residual power less than prev ue power
                UElist[ue_id].Set_Achievable_DataRate(0.0);
                curr += 1;
                break;
            }
            else { // residual power >= prev
                int mode = select_modulation_mode_base_on_residual_power(ap_id, ue_id, VLC_Max_Power-sum, 0, mode_table.size()-1, sum);

                if (mode == -1) { // condiction 2: cant find mode based on residual power
                    UElist[ue_id].Set_Achievable_DataRate(0.0);
                    continue;
                }
                else {
                    required_power = get_VLC_minimum_required_power(ap_id, ue_id, mode_table[mode], sum);
                    if (required_power < prev) { // condiction 3: required_power smaller than prev ue power
                        UElist[ue_id].Set_Achievable_DataRate(0.0);
                        continue;
                    }

                    // condiction 4: required_power safe
                    sum = sum + required_power;
                    UElist[ue_id].Set_Achievable_DataRate(mode_rate_table[mode]);
                    APlist[ap_id].Add_Associated_UE(ue_id, required_power);
                    curr += 1;
                    break;
                }
            }
        }
        else {
            UElist[ue_id].Set_Achievable_DataRate(mode_rate_table[UElist[ue_id].Get_Mode_Index()]);
            APlist[ap_id].Add_Associated_UE(ue_id, required_power);
            prev = required_power;
        }
    }

    for (; curr < tmp_UElist.size() ; curr++) {
        int ue_id = tmp_UElist[curr].GetID();
        UElist[ue_id].Set_Achievable_DataRate(0.0);
    }

    APlist[ap_id].Set_Residual_Power(VLC_Max_Power-sum);
}

void Update_RF_AP_power_allocation(int ap_id) {
    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();
    std::vector<UE_Node> tmp_UElist;
    for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
        tmp_UElist.push_back(UElist[UElist_of_AP[j]]);
    }

    APlist[ap_id].Clear_Associated_UE();
    sort(tmp_UElist.begin(), tmp_UElist.end(), compare_channel_gain);

    double sum = 0.0, prev = 0.0, required_power = 0.0;
    int curr = 0;

    for (; curr < tmp_UElist.size() ; curr++) {
        int ue_id = tmp_UElist[curr].GetID();

        required_power = get_RF_minimum_required_power(ap_id, ue_id, sum);
        required_power = (prev > required_power) ? prev : required_power;

        sum = sum + required_power;

        if (sum > RF_Max_Power) {
            sum = sum - required_power;
            if (RF_Max_Power - sum < prev) {
                UElist[ue_id].Set_Achievable_DataRate(0.0);
            }
            else {
                APlist[ap_id].Add_Associated_UE(ue_id, RF_Max_Power - sum);
                double SINR_based_on_residual_power = get_RF_SINR_based_on_Residual_Power(ap_id, ue_id, sum, RF_Max_Power-sum);
                UElist[ue_id].Set_Achievable_DataRate(get_RF_Achievable_Data_Rate(SINR_based_on_residual_power));
                sum = RF_Max_Power;
            }
            curr += 1;
            break;
        }
        else {
            UElist[ue_id].Set_Achievable_DataRate(UElist[ue_id].Get_Required_DataRate());
            APlist[ap_id].Add_Associated_UE(ue_id, required_power);
            prev = required_power;
        }
    }

    for (;curr < tmp_UElist.size() ; curr++) {
        int ue_id = tmp_UElist[curr].GetID();
        UElist[ue_id].Set_Achievable_DataRate(0.0);
    }

    APlist[ap_id].Set_Residual_Power(RF_Max_Power-sum);
}

void Update_All_AP_power_allocation() {
    for (int i = 0 ; i < APlist.size() ; i++) {
        if (i < RF_AP_Num)
            Update_RF_AP_power_allocation(i);
        else
            Update_VLC_AP_power_allocation(i);
    }
}

double generate_VHO_overhead() {
    std::normal_distribution<double> Gaussian = std::normal_distribution<double>(VHO_delay,0.05); // variance = 0.05s = 50ms
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
    return Gaussian(gen);
}

double generate_HHO_overhead() {
    std::normal_distribution<double> Gaussian = std::normal_distribution<double>(HHO_delay,0.05); // variance = 0.05s = 50ms
    std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
    return Gaussian(gen);
}


void Do_algorithm(){
    // 更新所有 UE 的位置與速度

    Update_UE_Condiction();

    /*UElist[0].SetPosition(-3.75, 3.75);
    UElist[0].Set_Required_DataRate(75.0);
    UElist[1].SetPosition(-3.6, 3.6);
    UElist[1].Set_Required_DataRate(75.0);
    UElist[2].SetPosition(-3.45, 3.45);
    UElist[2].Set_Required_DataRate(75);
    UElist[3].SetPosition(-3.3, 3.3);
    UElist[3].Set_Required_DataRate(75.0);
    Set_UE_default_modulation_mode();*/

    // 算所有 UE 的 Channel gain
    Calculate_Channel_Gain_Matrix(APlist, UElist, Channel_Gain_Matrix);

    // std::cout << bandwidth_per_cell * VLC_AWGN_spectral_density << std::endl;

    // std::cout << RF_AP_Bandwidth * RF_AWGN_spectral_density << std::endl;

    /*std::cout << UElist[0].Get_Required_DataRate()<< std::endl;
    std::vector<int> v1 = UElist[0].Get_Modulation_Mod();
    std::cout << v1[0] << v1[1] << v1[2] << v1[3] << std::endl;
    std::cout << UElist[1].Get_Required_DataRate()<< std::endl;
    std::vector<int> v2 = UElist[1].Get_Modulation_Mod();
    std::cout << v2[0] << v2[1] << v2[2] << v2[3] << std::endl;
    std::cout << UElist[2].Get_Required_DataRate()<< std::endl;
    std::vector<int> v3 = UElist[2].Get_Modulation_Mod();
    std::cout << v3[0] << v3[1] << v3[2] << v3[3] << std::endl;
    std::cout << UElist[3].Get_Required_DataRate()<< std::endl;
    std::vector<int> v4 = UElist[3].Get_Modulation_Mod();
    std::cout << v4[0] << v4[1] << v4[2] << v4[3] << std::endl;

    UElist[0].Set_Associated_AP(4);
    APlist[4].Add_Associated_UE(0, 0.0);
    UElist[1].Set_Associated_AP(4);
    APlist[4].Add_Associated_UE(1, 0.0);
    UElist[2].Set_Associated_AP(4);
    APlist[4].Add_Associated_UE(2, 0.0);
    //UElist[3].Set_Associated_AP(4);
    //APlist[4].Add_Associated_UE(3, 0.0);

    Update_VLC_AP_power_allocation(4);

    std::vector<int> AP4_UE_list = APlist[4].Get_UEs();
    for (int i = 0 ; i < AP4_UE_list.size() ; i++) {
        std::cout << "UE id " << AP4_UE_list[i] <<" " << Channel_Gain_Matrix[4][AP4_UE_list[i]] << " " << APlist[4].Get_UE_Power(AP4_UE_list[i]) << std::endl;
    }

    std::cout << "UE id 0 " << Channel_Gain_Matrix[4][0] << " " << UElist[0].Get_Achievable_DataRate() << std::endl;
    std::cout << "UE id 1 " << Channel_Gain_Matrix[4][1] << " " << UElist[1].Get_Achievable_DataRate() << std::endl;
    std::cout << "UE id 2 " << Channel_Gain_Matrix[4][2] << " " << UElist[2].Get_Achievable_DataRate() << std::endl;
    std::cout << "UE id 3 " << Channel_Gain_Matrix[4][3] << " " << UElist[3].Get_Achievable_DataRate() << std::endl;
    std::cout << "UE residual power " << APlist[4].Get_Residual_Power() << std::endl;
    std::cout << "UE residual power " << power_allocation_after_handover(4, 3)[1] << std::endl;*/

    if(Simulator::Now().GetSeconds() == 0) { // 初始化環境，UE 根據 channel gain 由大到小排序，採 SSS
        Set_UE_default_modulation_mode();
        Set_UE_Associated_AP_Based_On_Channel();

        std::vector<UE_Node> tmp_UElist;

        tmp_UElist.assign(UElist.begin(), UElist.end());

        sort(tmp_UElist.begin(), tmp_UElist.end(), compare_channel_gain);

        for (int i = 0 ; i < tmp_UElist.size() ; i++) {
            if (tmp_UElist[i].Get_Associated_AP() < RF_AP_Num) {

                // std::cout << tmp_UElist[i].Get_Associated_AP() << std::endl;
                // calculate required power
                int ap_id = tmp_UElist[i].Get_Associated_AP();
                int ue_id = tmp_UElist[i].GetID();
                double residual_power = APlist[ap_id].Get_Residual_Power();
                double prev = APlist[ap_id].Get_Prev_UE_Power();
                double required_power = get_RF_minimum_required_power(ap_id, ue_id, RF_Max_Power-residual_power);
                required_power = (prev > required_power) ? prev : required_power;

                std::cout << required_power << " " << residual_power << std::endl;
                // determine served by WiFi AP or not
                if (required_power <= residual_power) {
                    residual_power -= required_power;
                    APlist[ap_id].Set_Residual_Power(residual_power);
                    APlist[ap_id].Add_Associated_UE(ue_id, required_power);

                    UElist[ue_id].Set_Achievable_DataRate(UElist[ue_id].Get_Required_DataRate());
                } else {
                    UElist[ue_id].Set_Associated_AP(-1);
                    UElist[ue_id].Set_Achievable_DataRate(0.0);
                    UE_timer_mode[ue_id] = 1;
                }
            }
            else {
                // calculate required power
                int ap_id = tmp_UElist[i].Get_Associated_AP();
                int ue_id = tmp_UElist[i].GetID();
                double residual_power = APlist[ap_id].Get_Residual_Power();
                double prev = APlist[ap_id].Get_Prev_UE_Power();
                double required_power = get_VLC_minimum_required_power(ap_id, ue_id, VLC_Max_Power-residual_power);
                required_power = (prev > required_power) ? prev : required_power;

                // determine served by Lifi AP or not
                if (required_power <= residual_power) {
                    residual_power -= required_power;
                    APlist[ap_id].Set_Residual_Power(residual_power);
                    APlist[ap_id].Add_Associated_UE(ue_id, required_power);

                    UElist[ue_id].Set_Achievable_DataRate(mode_rate_table[UElist[ue_id].Get_Mode_Index()]);
                } else {
                    UElist[ue_id].Set_Associated_AP(-1);
                    UElist[ue_id].Set_Achievable_DataRate(0.0);
                    UE_timer_mode[ue_id] = 1;
                }
            }
        }
    }
    else {
        // update AP power
        Update_All_AP_power_allocation();

        for (int i = 0 ; i < UElist.size() ; i++) {
            if (UE_timer_mode[i] == 0) {
                if (UElist[i].GetVelocity() >= speed_threshold && UElist[i].Get_Associated_AP() >= RF_AP_Num) { // switch to RF AP immediately
                    UE_timer_mode[i] = 3;
                    UElist[i].Set_Achievable_DataRate(0.0);
                    Handover_timer[i].first = generate_VHO_overhead();
                    Handover_timer[i].second = Handover_timer[i].first;
                    Record_handover_time[i].first = Handover_timer[i].first;
                    Record_handover_time[i].second = Handover_timer[i].second;

                    int host_ap_id = UElist[i].Get_Associated_AP();
                    APlist[host_ap_id].Remove_Associated_UE(i);
                    Update_VLC_AP_power_allocation(host_ap_id);
                }
                else if (UElist[i].Get_Achievable_DataRate() == 0.0) {
                    UE_timer_mode[i] = 3;
                    Handover_timer[i].first = generate_VHO_overhead();
                    Handover_timer[i].second = std::min(Handover_timer[i].first, generate_HHO_overhead());
                    Record_handover_time[i].first = Handover_timer[i].first;
                    Record_handover_time[i].second = Handover_timer[i].second;
                }
                else if (UElist[i].Get_Achievable_DataRate() < UElist[i].Get_Required_DataRate()) {
                    UE_timer_mode[i] = 2;
                    Dwell_timer[i] = time_TTT;
                }
            }
            else if (UE_timer_mode[i] == 1) {
                int AP_num = 0;
                if (UElist[i].GetVelocity() >= speed_threshold)
                    AP_num = RF_AP_Num;
                else
                    AP_num = AP_Num;

                std::vector<std::pair<int, double>> AP_list;
                for (int j = 0 ; j < AP_num ; j++) {
                    std::vector<double> power = power_allocation_after_handover(j, i);
                    if (power[1] >= 0) {
                        AP_list.push_back({j, power[1]});
                    }
                }

                if (AP_list.size() == 0) { // dont have ap to handover, release ap resource
                    continue;
                } else {
                    sort(AP_list.begin(), AP_list.end(), compare_residual_power);

                    int ap_id = AP_list[0].first;
                    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();

                    int bigger_channel_ue = -1;
                    for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
                        if (Channel_Gain_Matrix[ap_id][UElist_of_AP[j]] > Channel_Gain_Matrix[ap_id][i]) {
                            bigger_channel_ue = j;
                            break;
                        }
                    }


                    if (bigger_channel_ue == -1) {
                        APlist[ap_id].Add_Associated_UE(i, 0.0);
                    } else {
                        APlist[ap_id].Insert_Associated_UE(i, bigger_channel_ue, 0.0);
                    }

                    UElist[i].Set_Associated_AP(ap_id);
                    if (ap_id < RF_AP_Num) {
                        Update_RF_AP_power_allocation(ap_id);
                    }
                    else {
                        Update_VLC_AP_power_allocation(ap_id);
                    }

                    UE_timer_mode[i] = 0;
                }
            }
            else if (UE_timer_mode[i] == 2) { // UE during dwell time
                if (UElist[i].GetVelocity() >= speed_threshold && UElist[i].Get_Associated_AP() >= RF_AP_Num) { // switch to RF AP immediately
                    UE_timer_mode[i] = 3;
                    Dwell_timer[i] = 0.0;
                    Handover_timer[i].first = generate_VHO_overhead();
                    Handover_timer[i].second = Handover_timer[i].first;
                    Record_handover_time[i].first = Handover_timer[i].first;
                    Record_handover_time[i].second = Handover_timer[i].second;

                    int host_ap_id = UElist[i].Get_Associated_AP();
                    APlist[host_ap_id].Remove_Associated_UE(i);

                    if (host_ap_id < RF_AP_Num)
                        Update_RF_AP_power_allocation(host_ap_id);
                    else
                        Update_VLC_AP_power_allocation(host_ap_id);
                }
                else if (Dwell_timer[i] > 0) { // during dwell time
                    if (UElist[i].Get_Achievable_DataRate() >= UElist[i].Get_Required_DataRate()) { // power recovery
                        UE_timer_mode[i] = 0;
                        Dwell_timer[i] = 0;
                    }
                    else if (UElist[i].Get_Achievable_DataRate() == 0) { // don't get any power
                        UE_timer_mode[i] = 3;
                        Dwell_timer[i] = 0;
                        Handover_timer[i].first = generate_VHO_overhead();
                        Handover_timer[i].second = std::min(Handover_timer[i].first, generate_HHO_overhead());
                        Record_handover_time[i].first = Handover_timer[i].first;
                        Record_handover_time[i].second = Handover_timer[i].second;
                    }
                    else { // still get some power, but cant get required power
                        Dwell_timer[i] -= delta_t;
                    }
                }
                else { // dwell time <= 0, if ue still can not serve by required data rate, doing handover.
                    if (UElist[i].Get_Achievable_DataRate() >= UElist[i].Get_Required_DataRate()) { // power recovery
                        UE_timer_mode[i] = 0;
                        Dwell_timer[i] = 0;
                    }
                    else {
                        UE_timer_mode[i] = 3;
                        Dwell_timer[i] = 0;
                        Handover_timer[i].first = generate_VHO_overhead();
                        Handover_timer[i].second = std::min(Handover_timer[i].first, generate_HHO_overhead());
                        Record_handover_time[i].first = Handover_timer[i].first;
                        Record_handover_time[i].second = Handover_timer[i].second;

                        int host_ap_id = UElist[i].Get_Associated_AP();
                        if (APlist[host_ap_id].Find_UE(i) == true) {// still get power from AP
                            APlist[host_ap_id].Remove_Associated_UE(i);

                            if (host_ap_id < RF_AP_Num)
                                Update_RF_AP_power_allocation(host_ap_id);
                            else
                                Update_VLC_AP_power_allocation(host_ap_id);
                        }
                    }
                }
            }
            else if (UE_timer_mode[i] == 3) {
                if (Handover_timer[i].first <= 0) { // can VHO
                    int host_ap_id = UElist[i].Get_Associated_AP();

                    int ap_id = 0, AP_num = 0;
                    if (UElist[i].GetVelocity() >= speed_threshold) {
                        AP_num = RF_AP_Num;
                    }
                    else
                        AP_num = AP_Num;

                    std::vector<std::pair<int, double>> AP_list;
                    for (; ap_id < AP_num ; ap_id++) {
                        std::vector<double> power = power_allocation_after_handover(ap_id, i);
                        if (power[1] >= 0) {
                            AP_list.push_back({ap_id, power[1]});
                        }
                    }

                    Handover_timer[i].first = 0.0;
                    Handover_timer[i].second = 0.0;

                    if (AP_list.size() == 0) { // dont have ap to handover, release ap resource
                        UElist[i].Set_Associated_AP(-1);
                        UE_timer_mode[i] = 1;
                        Record_handover_time[i].first = 0.0;
                        Record_handover_time[i].second = 0.0;
                    } else {
                        sort(AP_list.begin(), AP_list.end(), compare_residual_power);

                        ap_id = AP_list[0].first;
                        std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();

                        int bigger_channel_ue = -1;
                        for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
                            if (Channel_Gain_Matrix[ap_id][UElist_of_AP[j]] > Channel_Gain_Matrix[ap_id][i]) {
                                bigger_channel_ue = j;
                                break;
                            }
                        }

                        if (bigger_channel_ue == -1)
                            APlist[ap_id].Add_Associated_UE(i, 0.0);
                        else
                            APlist[ap_id].Insert_Associated_UE(i, bigger_channel_ue, 0.0);

                        UElist[i].Set_Associated_AP(ap_id);
                        if (ap_id < RF_AP_Num)
                            Update_RF_AP_power_allocation(ap_id);
                        else
                            Update_VLC_AP_power_allocation(ap_id);

                        UE_timer_mode[i] = 0;

                        if (host_ap_id == ap_id)
                            ;
                        else if (host_ap_id < RF_AP_Num && ap_id < RF_AP_Num) {
                            RFtoRF += 1;
                            HHO_count += 1;
                            total_HHO_delay += Record_handover_time[i].second;
                        }
                        else if (host_ap_id >= RF_AP_Num && ap_id < RF_AP_Num) {
                            VLCtoRF += 1;
                            VHO_count += 1;
                            total_VHO_delay += Record_handover_time[i].first;
                        }
                        else if (host_ap_id >= RF_AP_Num && ap_id >= RF_AP_Num) {
                            VLCtoVLC +=1;
                            HHO_count += 1;
                            total_HHO_delay += Record_handover_time[i].second;
                        }
                        else {
                            RFtoVLC +=1;
                            VHO_count += 1;
                            total_VHO_delay += Record_handover_time[i].first;
                        }

                        Record_handover_time[i].first = 0.0;
                        Record_handover_time[i].second = 0.0;
                    }

                }
                else if (Handover_timer[i].second <= 0) { // can HHO
                    int host_ap_id = UElist[i].Get_Associated_AP();

                    int ap_id = 0, AP_num = 0;
                    if (host_ap_id < RF_AP_Num) // host_ap is RF AP
                        ap_id = 0, AP_num = RF_AP_Num;
                    else
                        ap_id = RF_AP_Num, AP_num = AP_Num;

                    std::vector<std::pair<int, double>> AP_list;
                    for (; ap_id < AP_num ; ap_id++) {
                        std::vector<double> power = power_allocation_after_handover(ap_id, i);
                        if (power[1] >= 0) {
                            AP_list.push_back({ap_id, power[1]});
                        }
                    }

                    if (AP_list.size() == 0) { // dont have ap to handover, release ap resource
                        Handover_timer[i].first -= delta_t;
                        Handover_timer[i].second = 0.0;
                    } else {
                        sort(AP_list.begin(), AP_list.end(), compare_residual_power);

                        ap_id = AP_list[0].first;
                        std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();

                        int bigger_channel_ue = -1;
                        for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
                            if (Channel_Gain_Matrix[ap_id][UElist_of_AP[j]] > Channel_Gain_Matrix[ap_id][i]) {
                                bigger_channel_ue = j;
                                break;
                            }
                        }

                        if (bigger_channel_ue == -1) {
                            APlist[ap_id].Add_Associated_UE(i, 0.0);
                        } else {
                            APlist[ap_id].Insert_Associated_UE(i, bigger_channel_ue, 0.0);
                        }

                        UElist[i].Set_Associated_AP(ap_id);
                        if (ap_id < RF_AP_Num)
                            Update_RF_AP_power_allocation(ap_id);
                        else
                            Update_VLC_AP_power_allocation(ap_id);

                        Handover_timer[i].first = 0.0;
                        Handover_timer[i].second = 0.0;
                        UE_timer_mode[i] = 0;

                        if (host_ap_id == ap_id)
                            ;
                        else if (host_ap_id < RF_AP_Num && ap_id < RF_AP_Num)
                            RFtoRF += 1;
                        else if (host_ap_id >= RF_AP_Num && ap_id >= RF_AP_Num)
                            VLCtoVLC +=1;

                        HHO_count += 1;
                        total_HHO_delay += Record_handover_time[i].second;
                        Record_handover_time[i].first = 0.0;
                        Record_handover_time[i].second = 0.0;
                    }
                }
                else if (Handover_timer[i].first > 0 && Handover_timer[i].second > 0) { // still during handover
                    Handover_timer[i].first -= delta_t;
                    Handover_timer[i].second -= delta_t;
                }
            }

        }
    }

    update_analysis_data();

    std::cout << cnt << std::endl;
    std::cout << "============================" << std::endl;

    cnt += 1;
    Simulator::Schedule(Seconds(delta_t),&Do_algorithm);
}

void update_analysis_data() {
    for (int ue_id = 0 ; ue_id < UElist.size() ; ue_id++)
        UE_average_rate[ue_id] += UElist[ue_id].Get_Achievable_DataRate()/simulation_steps;

    int rf_ue_count = 0, vlc_ue_count = 0;
    for (int rf_ap_id = 0 ; rf_ap_id < RF_AP_Num ; rf_ap_id++)
        rf_ue_count += APlist[rf_ap_id].Associated_UE_Num();

    for (int vlc_ap_id = RF_AP_Num ; vlc_ap_id < AP_Num ; vlc_ap_id++)
        vlc_ue_count += APlist[vlc_ap_id].Associated_UE_Num();

    UE_connect_to_RF_AP += rf_ue_count;
    UE_connect_to_VLC_AP += vlc_ue_count;
    UE_connect_to_AP += (rf_ue_count + vlc_ue_count);
    Redundant_UE_count += (UE_Num - (rf_ue_count + vlc_ue_count));
}

void print_analysis_data() {
    double average_ue_number_connect_to_RF_AP = (double)UE_connect_to_RF_AP / simulation_steps;
    std::cout << "Average number of UE connect to RF AP: " << average_ue_number_connect_to_RF_AP << std::endl;

    double average_ue_number_connect_to_VLC_AP = (double)UE_connect_to_VLC_AP / simulation_steps;
    std::cout << "Average number of UE connect to VLC AP: " << average_ue_number_connect_to_VLC_AP << std::endl;

    double average_ue_number_connect_to_AP = (double)UE_connect_to_AP / simulation_steps;
    std::cout << "Average number of UE connect to AP: " << average_ue_number_connect_to_AP << std::endl;

    double average_number_redundant = (double)Redundant_UE_count / simulation_steps;
    std::cout << "Average number of redundant UE: " << average_number_redundant << std::endl;

    double average_ue_achievable_data_rate = 0.0, average_ue_required_data_rate = 0.0;
    for (int ue_id = 0 ; ue_id < UE_Num ; ue_id++) {
        average_ue_achievable_data_rate += UE_average_rate[ue_id];
        average_ue_required_data_rate += UElist[ue_id].Get_Required_DataRate();
    }

    average_ue_achievable_data_rate = average_ue_achievable_data_rate / UE_Num;
    average_ue_required_data_rate = average_ue_required_data_rate / UE_Num;
    std::cout << "Average achievable data rate of UE: " << average_ue_achievable_data_rate << std::endl;
    std::cout << "Average required data rate of UE: " << average_ue_required_data_rate << std::endl;

    //=======handover=======

    std::cout << "Number of VHO: " << VHO_count << std::endl;
    std::cout << "Number of RF to VLC: " << RFtoVLC << std::endl;
    std::cout << "Number of VLC to RF: " << VLCtoRF << std::endl;
    std::cout << "Number of HHO: " << HHO_count << std::endl;
    std::cout << "Number of VLC to VLC: " << VLCtoVLC << std::endl;
    std::cout << "Number of RF to RF: " << RFtoRF << std::endl;
}
