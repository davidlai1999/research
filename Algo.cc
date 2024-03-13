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

std::vector<double> BER_required_SINR(number_of_layer+1,0.0);

//channel gain
std::vector<std::vector<double>> Channel_Gain_Matrix(RF_AP_Num+VLC_AP_Num,std::vector<double> (UE_Num,0));

//frequency reuse
std::vector<int> frequency (VLC_AP_Num,0);

//dwell time timer
std::vector<double> Dwell_timer (UE_Num, 0.0);

//timer of handover, first is VHO, second is HH0
std::vector<std::pair<double, double>> Handover_timer (UE_Num, {0.0, 0.0});

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
std::vector<int> Redundant_UElist;

struct less_than_key
{
    inline bool operator() (const UE_Node& node1, const UE_Node& node2)
    {
        return (Channel_Gain_Matrix[node1.Get_Associated_AP()][node1.GetID()] > Channel_Gain_Matrix[node2.Get_Associated_AP()][node2.GetID()]);
    }
};

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
        Ptr<MobilityModel> AP_MobilityModel = (AP_Nodes.Get(i))->GetObject<MobilityModel> ();
        Vector pos = AP_MobilityModel->GetPosition ();
        APlist.push_back(AP_Node(i, pos, AP_Nodes.Get(i)));
    }

    //frequency reuse setting
    frequency[1]=frequency[7]=frequency[9]=frequency[15]=1;
    frequency[2]=frequency[4]=frequency[10]=frequency[12]=2;
    frequency[3]=frequency[5]=frequency[11]=frequency[13]=3;
}

void Update_UE_Condiction(std::vector<UE_Node> & UElist) {
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
    return bandwidth_per_rb / pow(2, layer) * log2(1 + SINR);
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
    double AWGN = bandwidth_per_rb * AWGN_Power_spectral_density;
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

}

double get_RF_minimum_required_power(const int ap_id, const int ue_id, const double sum_prev_ue_power) {
    double path_loss = 0.0;

    double d = GetDistance_AP_UE(APlist[ap_id].Get_Node(),UElist[ue_id].Get_Node());

    if(d < d_ref)
        path_loss = 20 * log10(central_freq * d) - 147.5;
    else
        path_loss = 20 * log10(central_freq * pow(d,2.75)/pow(d_ref,1.75)) - 147.5;

    double channel_gain = path_loss * std::pow(Channel_Gain_Matrix[ap_id][ue_id], 2);

    double result = channel_gain * sum_prev_ue_power;

    result = result + bandwidth_per_rb * AWGN_Power_spectral_density;

    result = result * (pow(2, UElist[ue_id].Get_Required_DataRate() / bandwidth_per_rb) - 1);

    result = result / channel_gain;
}

std::vector<double> power_allocation_after_handover(int ap_id, int ue_id) {
    if (APlist[ap_id].Get_Residual_Power() <= 0) return {0.0, APlist[ap_id].Get_Residual_Power()};
    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();

    double sum = 0.0, prev = 0.0, required_rate = 0.0;
    int bigger_channel_ue = -1;
    for (int i = 0 ; i < UElist_of_AP.size() ; i++) {
        if (Channel_Gain_Matrix[ap_id][UElist_of_AP[i]] <= Channel_Gain_Matrix[ap_id][ue_id]) {
            prev = APlist[ap_id].Get_UE_Power(UElist_of_AP[i]);
            sum = sum + prev;
        }
        else {
            bigger_channel_ue = i;
            break;
        }
    }

    if (ap_id < RF_AP_Num)
        required_rate = get_RF_minimum_required_power(ap_id, ue_id, sum);
    else
        required_rate = get_VLC_minimum_required_power(ap_id, ue_id, sum);

    required_rate = (prev > required_rate) ? prev : required_rate;
    sum = sum + required_rate;
    prev = required_rate;

    if (bigger_channel_ue != -1) {
        for (int i = bigger_channel_ue ; i < UElist_of_AP.size() ; i++) {
            if (ap_id < RF_AP_Num) {
                required_rate = get_RF_minimum_required_power(ap_id, ue_id, sum);
                required_rate = (prev > required_rate) ? prev : required_rate;
                sum = sum + required_rate;
                prev = required_rate;
            }
            else {
                required_rate = get_VLC_minimum_required_power(ap_id, ue_id, sum);
                required_rate = (prev > required_rate) ? prev : required_rate;
                sum = sum + required_rate;
                prev = required_rate;
            }

        }
    }

    double residual_power = VLC_Max_Power - sum;
    return {required_rate, residual_power};
}

int find_Best_VLC_Channel_Gain(int ue_id) {
    int result = -1;
    double maxi = 1e-9;
    for (int ap_id = RF_AP_Num ; ap_id < VLC_AP_Num ; ap_id++) {
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

void set_UE_Associated_AP_Based_On_Channel() {
    for (int i = 0 ; i < UElist.size() ; i++) {
        int ap_id = find_Best_VLC_Channel_Gain(i);

        if (ap_id == -1)
            ap_id = find_Best_RF_Channel_Gain(i);

        UElist[i].Set_Associated_AP(ap_id);
    }
}

void Update_AP_power_allocation(int ap_id) {
    std::vector<int> UElist_of_AP = APlist[ap_id].Get_UEs();
    std::vector<UE_Node> tmp_UElist;
    for (int j = 0 ; j < UElist_of_AP.size() ; j++) {
        tmp_UElist.push_back(UElist[UElist_of_AP[j]]);
    }

    APlist[ap_id].Clear_Associated_UE();
    sort(tmp_UElist.begin(), tmp_UElist.end(), less_than_key());

    double sum = 0.0, prev = 0.0;
    double required_power = 0.0;
    for (int j = 0 ; j < tmp_UElist.size() ; j++) {
        int ue_id = tmp_UElist[j].GetID();

        if (ap_id < RF_AP_Num)
            required_power = get_RF_minimum_required_power(ap_id, ue_id, sum);
        else // is VLC AP
            required_power = get_VLC_minimum_required_power(ap_id, ue_id, sum);

        required_power = (prev > required_power) ? prev : required_power;
        APlist[ap_id].Add_Associated_UE(ue_id, required_power);

        sum = sum + required_power;
        prev = required_power;
    }

    APlist[ap_id].Set_Residual_Power(VLC_Max_Power-sum);
}

void Update_All_AP_power_allocation() {
    for (int i = 0 ; i < APlist.size() ; i++) {
        Update_AP_power_allocation(i);
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
    Update_UE_Condiction(UElist);

    /*UElist[0].SetPosition(-3.75, 3.75);
    UElist[1].SetPosition(-3.2, 3.2);
    UElist[2].SetPosition(-3.0, 3.0);
    UElist[3].SetPosition(-2.7, 2.7);*/

    // 算所有 UE 的 Channel gain
    Calculate_Channel_Gain_Matrix(APlist, UElist, Channel_Gain_Matrix);

    /*int mode = select_modulation_mode(UElist[0].Get_Required_DataRate(), 0, mode_rate_table.size());
    UElist[0].Set_Modulation_Mod(mode_table[mode]);
    mode = select_modulation_mode(UElist[1].Get_Required_DataRate(), 0, mode_rate_table.size());
    UElist[1].Set_Modulation_Mod(mode_table[mode]);
    mode = select_modulation_mode(UElist[2].Get_Required_DataRate(), 0, mode_rate_table.size());
    UElist[2].Set_Modulation_Mod(mode_table[mode]);
    mode = select_modulation_mode(UElist[3].Get_Required_DataRate(), 0, mode_rate_table.size());
    UElist[3].Set_Modulation_Mod(mode_table[mode]);
    std::cout << UElist[0].Get_Required_DataRate()<< std::endl;
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

    // std::cout << APlist[4].GetPosition().x << APlist[4].GetPosition().y << std::endl;

    double prev = get_VLC_minimum_required_power(4, 0, 0.0), sum = 0.0;
    std::cout << Channel_Gain_Matrix[4][0] << " " << prev << std::endl;
    //std::cout << get_channel(APlist[4], UElist[0]) << std::endl;
    sum = sum + prev;
    prev = get_VLC_minimum_required_power(4, 1, sum);
    std::cout << Channel_Gain_Matrix[4][1] << " " << prev << std::endl;
    sum = sum + prev;
    prev = get_VLC_minimum_required_power(4, 2, sum);
    std::cout << Channel_Gain_Matrix[4][2] << " " << prev << std::endl;
    sum = sum + prev;
    prev = get_VLC_minimum_required_power(4, 3, sum);
    std::cout << Channel_Gain_Matrix[4][3] << " " << prev << std::endl;
    sum = sum + prev;
    std::cout << sum << std::endl;
    if (sum >= VLC_Max_Power) {
        sum = sum - prev;
        std::cout << sum << std::endl;
        int mode = select_modulation_mode_base_on_residual_power(4, 3, VLC_Max_Power-sum, 0, mode_table.size()-1, sum);

        prev = get_VLC_minimum_required_power(4, 3, mode_table[mode], sum);
        std::cout << VLC_Max_Power-sum << std::endl;
        std::cout << mode << " " << prev << std::endl;
    }
    std::cout << "========================" << std::endl;*/

    if(Simulator::Now().GetSeconds() == 0) { // 初始化環境，UE 根據 channel gain 由大到小排序，採 SSS
        set_UE_Associated_AP_Based_On_Channel();

        std::vector<UE_Node> tmp_UElist;

        tmp_UElist.assign(UElist.begin(), UElist.end());

        sort(tmp_UElist.begin(), tmp_UElist.end(), less_than_key());

        for (int i = 0 ; i < tmp_UElist.size() ; i++) {
            if (tmp_UElist[i].Get_Associated_AP() < RF_AP_Num) {
                // std::cout << tmp_UElist[i].Get_Associated_AP() << std::endl;
                // calculate required power
                int ue_id = tmp_UElist[i].GetID();
                int ap_id = find_Best_RF_Channel_Gain(ue_id);
                double residual_power = APlist[ap_id].Get_Residual_Power();
                double required_power = get_RF_minimum_required_power(ap_id, ue_id, VLC_Max_Power-residual_power);

                // determine served by WiFi AP or not
                if (required_power <= residual_power) {
                    residual_power -= required_power;
                    APlist[ap_id].Set_Residual_Power(residual_power);
                    APlist[ap_id].Add_Associated_UE(ue_id, required_power);
                } else {
                    UElist[ue_id].Set_Associated_AP(-1);
                    Redundant_UElist.push_back(ue_id);
                }
            }
            else {
                // select modulation mode
                int mode = select_modulation_mode(tmp_UElist[i].Get_Required_DataRate(), 0, mode_rate_table.size()-1);
                tmp_UElist[i].Set_Modulation_Mod(mode_table[mode]);

                // calculate required power
                int ap_id = tmp_UElist[i].Get_Associated_AP();
                int ue_id = tmp_UElist[i].GetID();
                double residual_power = APlist[ap_id].Get_Residual_Power();
                double required_power = get_RF_minimum_required_power(ap_id, ue_id, VLC_Max_Power-residual_power);

                if (APlist[ap_id].Associated_UE_Num() > 0) {
                    // check if pre ue required power is bigger than current ue
                }

                // determine served by Lifi AP or not
                if (required_power <= residual_power) {
                    residual_power -= required_power;
                    APlist[ap_id].Set_Residual_Power(residual_power);
                    APlist[ap_id].Add_Associated_UE(ue_id, required_power);
                } else {
                    UElist[ue_id].Set_Associated_AP(-1);
                    Redundant_UElist.push_back(ue_id);
                }
            }
        }
    }
    else {
        // update AP power
        Update_All_AP_power_allocation();

        for (int i = 0 ; i < UElist.size() ; i++) {
            if (UE_timer_mode[i] == 0) {
                ;
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

                    Update_AP_power_allocation(ap_id);
                    UElist[i].Set_Associated_AP(ap_id);
                }
            }
            else if (UE_timer_mode[i] == 2) { // UE during dwell time
                if (UElist[i].GetVelocity() >= speed_threshold) { // switch to RF AP immediately
                    UE_timer_mode[i] = 3;
                    Dwell_timer[i] = 0.0;
                    Handover_timer[i].first = generate_VHO_overhead();
                    Handover_timer[i].second = std::min(Handover_timer[i].first, generate_HHO_overhead());

                    int host_ap_id = UElist[i].Get_Associated_AP();
                    APlist[host_ap_id].Remove_Associated_UE(i);
                    Update_AP_power_allocation(host_ap_id);
                }
                else if (Dwell_timer[i] > 0) {
                    // test ue speed
                    ;
                }
                else { // dwell time <= 0, if ue still can not serve by required data rate, doing handover.
                    ;
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

                        Update_AP_power_allocation(ap_id);
                        UElist[i].Set_Associated_AP(ap_id);
                        UE_timer_mode[i] = 0;
                    }

                }
                else if (Handover_timer[i].second <= 0) { // can HHO
                    int host_ap_id = UElist[i].Get_Associated_AP();

                    int ap_id = 0, AP_num = 0;
                    if (host_ap_id < RF_AP_Num) {// host_ap is RF AP
                        ap_id = 0, AP_num = RF_AP_Num;
                    }
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

                        Update_AP_power_allocation(ap_id);
                        UElist[i].Set_Associated_AP(ap_id);
                        Handover_timer[i].first = 0.0;
                        Handover_timer[i].second = 0.0;
                        UE_timer_mode[i] = 0;
                    }
                }
                else if (Handover_timer[i].first > 0 && Handover_timer[i].second > 0) { // still during handover
                    Handover_timer[i].first -= delta_t;
                    Handover_timer[i].second -= delta_t;
                }
            }

        }
    }

    /*for (int i = 0 ; i < APlist.size() ; i++) {
        std::cout << APlist[i].Get_Residual_Power() << std::endl;
    }*/

    Simulator::Schedule(Seconds(delta_t),&Do_algorithm);
}

