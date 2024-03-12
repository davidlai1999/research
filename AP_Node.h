#ifndef AP_NODE_H
#define AP_NODE_H

#include <iostream>
#include <fstream>
#include <string>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "global_environment.h"

using namespace ns3;

class AP_Node {

    public :
        AP_Node(int id,Vector pos, Ptr<Node> node); //constructor

        void Set_AP_ID(int id); //set ap id

        int GetID(void); // get ap id

        void SetPosition(Vector pos_from_model); // set ap position

        Vector GetPosition(void); // get ap position

        void Add_Associated_UE(int uid, double power); // add an UE to the AP

        void Insert_Associated_UE(int uid, int dst, double power);

        void Remove_Associated_UE(int uid); // remove an UE from the AP

        void Clear_Associated_UE(); // remove all UE from the AP

        std::vector<int> Get_UEs(void); // get all the UEs under an AP

        Ptr<Node> Get_Node(void); //get the ns3 node of the AP

        bool Find_UE(int uid); // is the UE under this AP?

        void Set_UE_Power(int uid, double required_power); // set the required resource needed to satisfy the UE's required rate and the actual resource that is assigned to the UE

        double Get_UE_Power(int uid); // get the actual power of an UE

        double Get_Required_Power(int uid); // get the required power of an UE

        double Get_Residual_Power(void); // get the residual power of the AP

        void Set_Residual_Power(double Power); // set the residual power of the AP
    private :
        int Node_ID;                                             //ap id
        Vector pos;                                              //ap position
        Ptr<Node> node;                                          //ns3 node
        double residual_power;                                   //residual power
        std::vector<std::pair<int, double>> accociated_UEs; //list of the associated UEs and their power
};
#endif

