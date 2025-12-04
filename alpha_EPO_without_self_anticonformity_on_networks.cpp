#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <iterator>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
using namespace std;

struct Node {
    int pub;
    int priv;
    int num_neighbors;
    vector <int> neighbors;
};


double c_priv (Node* network, int N) {
    double c = 0;
    for (int i = 0; i < N; i++){
        c += (1 + network[i].priv);
    }
    return c/2/(double)N;
}

double c_pub (Node* network, int N) {
    double c = 0;
    for (int i = 0; i < N; i++){
        c += (1 + network[i].pub);
    }
    return c/2/(double)N;    
}

double dissonance (Node* network, int N) {
    double d = 0;
    for (int i = 0; i < N; i++){
        d += (1 - network[i].pub*network[i].priv);
    }
    return d/(double)(2*N);   
}

void exp_results (int N, int q, double alpha, string ntw_name,  double p, double pub, double priv, double diss){
  string fileName = "alpha_EPO_withoutSA_N" + to_string(N)+ "_q" + to_string(q)+"_alpha"+to_string((int)(alpha*100))+"_"+ntw_name+".txt";

  ofstream resFile;
  resFile.open(fileName, ios_base::app);
  resFile<<fixed<<setprecision(3)<<p<<"\t"<<fixed<<setprecision(3)<<pub<<"\t"<<fixed<<setprecision(3)<<priv<<"\t"<<fixed<<setprecision(3)<<diss<<endl;
  resFile.close();
};


int main()
{
    auto start = chrono::high_resolution_clock::now();
    int N = 0;
    int q = 3;
    int q_temp;
    int MCS = 5000;
    int MCS_term = 3000;
    int step = 1;
    int number_of_networks = 1;
    double f = 0.5;
    // double alpha = 0.1;
    double alpha_min = 0.1;
    double alpha_max = 0.9; 
    double alpha_step = 0.1;
    double p_max = 1.;
    double p_step = 0.01;
    vector <int> nodes;
    vector <int> q_voters;   
    vector <int> neighbors_;
    string network_name;
    vector <string> networks_names = {"random_graph_N_10000_k_15", "random_graph_N_10000_k_50", "random_graph_N_10000_k_150"};
    // {"WS_N_10000_k_15_c_0.10", "WS_N_10000_k_50_c_0.10", "WS_N_10000_k_150_c_0.10",   "WS_N_10000_k_150_c_0.30", "WS_N_10000_k_150_c_0.50" };
    //{"WS_N_100000_k_15_c_0.10", "WS_N_100000_k_50_c_0.10", "WS_N_100000_k_150_c_0.10", "WS_N_100000_k_15_c_0.30", "WS_N_100000_k_50_c_0.30", "WS_N_100000_k_150_c_0.30", "WS_N_100000_k_15_c_0.50", "WS_N_100000_k_50_c_0.50", "WS_N_100000_k_150_c_0.50" };
    //{"WS_N_10000_k_15_c_0.50", "WS_N_10000_k_15_c_0.30", "WS_N_10000_k_50_c_0.30", "WS_N_10000_k_50_c_0.50"}; 
    //{"S1", "S2", "M1", "M2", "L1", "L2"};
    //   

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    auto gen = std::mt19937(seed);
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);
    uniform_int_distribution<int> int_distribution(0, 1);
    random_device rd;
    mt19937 g(rd());

    for (int nn = 0; nn < networks_names.size(); nn++) {
        //network_name = "WS_as_" + networks_names[nn] ;
        network_name = networks_names[nn];
        string network_file =  "..\\Networks_random_graph\\"+network_name + ".txt";
        // "..\\Networks_WS\\"+network_name + ".txt";
        //string network_file = "..\\Networks_from_organizations\\"+network_name + ".txt";
        N = 0;
        ifstream file(network_file);
        // String to store each line of the file.
        string line;        
        if (file.is_open()) {
                // Read each line from the file and store it in the
                // 'line' variable.
                while (getline(file, line)) {
                    N++;
                }

                // Close the file stream once all lines have been
                // read.
                file.close();
                cout<<"no. nodes = "<<N<<endl;
            }
            else {
                // Print an error message to the standard error
                // stream if the file cannot be opened.
                cerr << "Unable to open file!" << endl;
            }
        uniform_int_distribution<int> nodes_distribution(0, N - 1);

        for (double alpha = alpha_min; alpha <= alpha_max +0.0001; alpha += alpha_step) {
            cout<<network_name<<'\t'<<"alpha = "<<alpha<<'\n'<<endl;
            for (double p = 0; p <= p_max + 0.001; p += p_step){
                double c_p_priv = 0;
                double c_p_pub = 0;
                double diss_p = 0; 
                double d_A_p = 0;
                double d_T_p = 0;
                double number_of_mes = 0;
                double number_of_mes_A = 0;        
                double number_of_mes_T = 0;
                for (int z = 0; z < number_of_networks; z++) {
                    Node* network = new Node[N];  


                ifstream file(network_file);
                // String to store each line of the file.
                string line;
                if (file.is_open()) {
                    // Read each line from the file and store it in the
                    // 'line' variable.
                    int idx = 0;
                    while (getline(file, line)) {
                        std::stringstream iss(line);
                        std::vector<int> neighbors_;
                        int v;
                        //cout<<ss<<endl;
                        while (iss >> v)                 // populate the new vector
                        {
                            neighbors_.push_back(v);
                            //cout<<v<<endl;
                        }
                        network[idx] = Node{1, 1, (int)neighbors_.size(), neighbors_};
                        //cout<<(int)neighbors_.size()<<endl;
                        idx++;
                    }

                    // Close the file stream once all lines have been
                    // read.
                    file.close();
                }
                else {
                    // Print an error message to the standard error
                    // stream if the file cannot be opened.
                    cerr << "Unable to open file!" << endl;
                }



                    //cout<<"network done"<<endl;

                    // Monte Carlo steps            
                    for (int mcs = 0; mcs < MCS; mcs++) {
                        //cout<<mcs<<endl;
                        int node_index;
                        double r;

                        //single step
                        for (int i = 0; i < N; i++) {
                            node_index = nodes_distribution(generator);
                            if (distribution(generator) > alpha){
                                //public update
                                
                                r = distribution(generator);
                                if (r < p) {
                                    network[node_index].pub = network[node_index].priv;
                                } else {
                                    q_temp = min(network[node_index].num_neighbors, q);
                                    sample(network[node_index].neighbors.begin(), network[node_index].neighbors.end(), back_inserter(q_voters), q_temp, gen);
                                    
                                    //disinhibitory contagion
                                    if (network[node_index].priv != network[node_index].pub){
                                        for (int i = 0; i < q_temp; i++){
                                            if (network[q_voters[i]].pub == network[node_index].priv){
                                                network[node_index].pub = network[node_index].priv;
                                                break;
                                            }
                                        }
                                    } else {
                                    //compliance
                                        int opinion = network[q_voters[0]].pub;
                                        bool unamity = true; 
                                        for (int i = 0; i < q_temp; i++){
                                            if (network[q_voters[i]].pub != opinion){
                                                unamity = false; 
                                                break;
                                            } 
                                        }
                                        if (unamity){
                                            network[node_index].pub = opinion;
                                        }
                                    }      

                                    q_voters.clear();      
                                    
                                }

                            } else {
                                
                                //private update
                                
                                r = distribution(generator);
                                if (r < p) {
                                    if (distribution(generator) < f ){network[node_index].priv *= (-1);}                     
                                } else {
        
                                    q_temp = min(network[node_index].num_neighbors, q);                               
                                    sample(network[node_index].neighbors.begin(), network[node_index].neighbors.end(), back_inserter(q_voters), q_temp, gen);
                                    int opinion = network[node_index].pub;
                                    bool unamity = true; 
                                        for (int i = 0; i < q_temp; i++){
                                            if (network[q_voters[i]].pub != opinion){
                                                unamity = false; 
                                                break;
                                            } 
                                        }
                                        if (unamity){
                                            network[node_index].priv = opinion;
                                        }
                                    
                                    q_voters.clear();
                                
                                                
                                }                    
                            
                            }

                            
                        }

                        if (mcs % step == 0 && mcs > MCS_term) {
                            double c_priv_ = c_priv(network,  N);
                            double c_pub_ = c_pub(network,  N);
                            double d_ = dissonance(network,  N);
                            c_p_priv += c_priv_;
                            c_p_pub += c_pub_;
                            diss_p += d_;
                            number_of_mes++;
                        }

                        

                    }

                    
                    
                    delete[] network;
                }
                c_p_priv = c_p_priv / number_of_mes;
                c_p_pub = c_p_pub / number_of_mes;
                diss_p = diss_p / number_of_mes;

                std::cout<<p<<"\t"<<c_p_pub<<"\t"<<c_p_priv<<"\t"<<diss_p<<"\t"<<d_A_p<<"\t"<<d_T_p<< endl;
                exp_results(N, q, alpha,network_name, p, c_p_pub, c_p_priv, diss_p);
            }    
            
                
            auto finish = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = finish - start;
            std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        }
    }

}
