#ifndef BRAIN_HPP
#define BRAIN_HPP
#include <iostream>
#include <random>
#include <vector>
#include <deque>
#include <set>
#include <unordered_set>
#include <string>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <cmath>
#include <boost/math/distributions/binomial.hpp>
#include <chrono>

using std::cout;
using std::deque;
using std::endl;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

#include "truncated_normal.cpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <fstream>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/container/flat_set.hpp> //需要能快速查找，并快速随机取样所以用这个

using boost::container::flat_set;
using boost::math::policies::discrete_quantile;
using boost::math::policies::integer_round_up;
using boost::math::policies::policy;
using boost::math::policies::real;

// typedef boost::math::binomial_distribution<double, policy<discrete_quantile<real>>> binom_real_quantile;
typedef boost::math::binomial_distribution<double, policy<discrete_quantile<integer_round_up>>> binom_round_up;

class Stimulus
{
public:
    int k; // A stimulus is a set of about k neurons firing simultaneously (“presented”) in the sensory area.
    Stimulus(int k0 = 0);
};

class Area
{
public:
    string name; // Area ID
    //在python的代码中，name应该是一个字符串，通过dict查找area对象

    int n;                                      //该脑区中神经元的个数 Number of neurons
    int k;                                      //一轮中被激活的神经元个数 Number of winners
    float beta;                                 // Default plasticity
    unordered_map<string, float> stimulus_beta; // Plasticity for connections from stimuli to this area
    unordered_map<string, float> area_beta;     // Plasiticity for connections from areas to this area
    int w;                                      // width, Number of winners, w是该区域上一轮的winner对这一轮输入的数量,也是这一轮中有连接记录的神经元数量
    flat_set<int> winners;                      // List of winners
    int new_w;                                  // Number of new winners,
    flat_set<int> new_winners;                  // List of new winners, 保存了所有在本轮被激活神经元的编号
    deque<int> saved_winners;                   // List of lists of winners at each round
    vector<int> saved_w;                        // List of numbers of winners at each round 用于记录历史的w
    int num_first_winners;                      // Initial number of winners init -1,这轮结束后统计的新出现记录的神经元数量
    bool fixed_assembly;                        // Whether to fix (freeze) the assembly (winners) in this area init false
    bool explicit_;                             // Whether to run the simulation explicitly init false

    Area(string name0 = "NULL", int n0 = 0, int k0 = 0, float beta0 = 0.05);
    void update_winners();
    void update_stimulus_beta(string name0, float new_beta0);
    void update_area_beta(string name0, float new_beta0);
    void fix_assembly();
    void unfix_assembly();
};

struct compare_col_first
{
    bool operator()(const pair<int, int> &p1, const pair<int, int> &p2) const
    {
        bool result = false;
        if (p1.second < p2.second)
        {
            result = true;
        }
        else if (p1.second == p2.second && p1.first < p2.first)
        {
            result = true;
        }
        return result;
    }
};

class Matrix
{ // to save connectome data
private:
    map<pair<int, int>, float, compare_col_first> data; //用元组的方式保存非零元素
    // default: sort by pair.first then pair.second, from small to large

    int max_row;
    int max_col;

public:
    Matrix();
    float &at(int row, int col) { return this->data[pair<int, int>(row, col)]; }         // do not check index, if not exist, return 0;
    float &real_at(int row, int col) { return this->data.at(pair<int, int>(row, col)); } // do check index
    float read_at(int row, int col);
    void multiply_at(int row, int col,const float& rate);//write if there is element
    float get_col_winner_sum(int col, flat_set<int> &winners); // for the col, sum row in winners
    int get_max_row() { return max_row; }
    int get_max_col() { return max_col; }
    void set_row(int new_row_size);
    void set_col(int new_col_size);
    void dump(string file_name); // dump and free original
    void load(string file_name);
    map<pair<int, int>, float>::iterator begin() { return data.begin(); }
    map<pair<int, int>, float>::iterator end() { return data.end(); }
    void print();
    void clear();
    void set_shape(int row, int col);
};

class Brain
{
public:
    unordered_map<string, Area> areas;                                              // Dictionary of areas in the brain                                      //考虑到脑区编号无意义，创建时从0开始从小到大使用，可用deque代替unordered_map/dict
    unordered_map<string, Stimulus> stimuli;                                        // Dictionary of stimuli in the brain
    unordered_map<string, unordered_map<string, deque<float>>> stimuli_connectomes; // For each stimulus in the brain, dictionary of connections from it to each area of the brain
    unordered_map<string, unordered_map<string, class Matrix>> connectomes;
    float p; // probability parameter for the edges of the underlying Erdos-Renyi graph
    bool save_size;
    bool save_winners;
    bool no_plasticity; // init false # For debugging purposes in applications (eg. language)

    Brain(float p0 = 0.01, bool save_size0 = true, bool save_winners0 = false);

    // # Add stimulus to brain. The weights of the connections to an explicit area of the brain are set randomly
    // # according to a binomial distribution with parameters k and the value p of the brain. The plasticity of
    // # these connections are set equal to the default plasticity of the area.
    void add_stimulus(string name0, int k0);

    // # Add a non-explicit area to the brain. Since the area is not explicit, all connections from the
    // # stimuli of the brain and from/to all the areas of the brain are initially set to empty (they
    // # will be set during the project phase in order to improve performance).
    void add_area(string name0, int n0, int k0, float beta0);

    // # Add an explicit area to the brain. Since the area is explicit, the weights of all connections
    // # from a stimuli of the brain the new area are initially set randomly according to a
    // # binomial distribution with parameters the value k of the stimulus and the value p of the brain
    // # (with the default plasticity). The weights of all connections from/to an explicit area of the brain
    // # to the new area are initially and fully set randomly according to a binomial distribution
    // # with parameters 1 and the value p of the brain. The weights of all connections from/to a
    // # non-explicit area of the brain to the new area are initially set to empty.
    // # In all cases, the plasticity of the connections is set to the default plasticity.
    // # The number of winners of the new area is set equal to the number of its neurons.
    void add_explicit_area(string name0, int n0, int k0, int beta0);

    // # Update the plasticities of the connections between stimuli and areas. Each area update consists of
    // # of the destination area and a list of update rules: each update rule specifies the source area
    // # and the new plasticity. Each stimulus update consists of the destination area and a list of update rules:
    // # each update rule specifies the source stimulus and the new plasticity.
    // 可以考虑使用unordered_multiunordered_map
    void update_plasticities(unordered_map<string, unordered_map<string, float>> area_update_unordered_map, unordered_map<string, unordered_map<string, float>> stim_update_unordered_map);
    // area_update_unordered_map consists of area1: list[ (area2, new_beta) ]
    // represents new plasticity FROM area2 INTO area1
    // stim_update_unordered_map consists of area: list[ (stim, new_beta) ]
    // represents new plasticity FROM stim INTO area

    // # Execute the project from stimuli and/or areas to areas. For each stimulus (key) in the first dictionary,
    // # the list (value) of areas to which the stimulus has the project is specified. For each area (key),
    // # in the second dictionary, the list (value) of areas to which the area has the project is specified.
    // # The function collects, for each area, the list of stimuli and areas that project to it (basically, it
    // # computes the inverse of the input mappings). Then, for each area which has to "receive" a projection
    // # (from either stimuli or areas), it invokes the function which actually performs the projection (this
    // # function returns the number of winners in the destination area). If the new winners have to be later
    // # analysed, then their list is appended to the the list of lists of winners of the area. When everything
    // # has been done, the function updates the destination areas.
    void project(const unordered_map<string, deque<string>> &stim_to_area, const unordered_map<string, deque<string>> &area_to_area, bool verbose = false);
    // Validate stim_area, area_area well defined
    // stim_to_area: {"stim1":["A"], "stim2":["C","A"]}
    // area_to_area: {"A":["A","B"],"C":["C","A"]}

    int project_into(Area &area0, deque<string> from_stimuli, deque<string> from_areas, bool verbose = false);
    // projecting everything in from stim_in[area] and area_in[area]
    // calculate: inputs to self.connectomes[area] (previous winners)
    // calculate: potential new winners, Binomial(sum of in sizes, k-top)
    // k top of previous winners and potential new winners
    // if new winners > 0, redo connectome and intra_connectomes
    // have to wait to replace new_winners
};

#include "brain.cpp"
#endif