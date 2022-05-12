#ifndef BRAIN_CPP
#define BRAIN_CPP

#include "brain.hpp"

using Stimulus = class Stimulus;
using Area = class Area;
using Brain = class Brain;
using Matrix = class Matrix;

std::random_device rd{};
std::mt19937 generator{rd()};

Matrix::Matrix()
{
    this->max_row = 0;
    this->max_col = 0;
}

void Matrix::set_row(int new_row_size)
{
    this->max_row = new_row_size;
}
void Matrix::set_col(int new_col_size)
{
    this->max_col = new_col_size;
}

void Matrix::dump(string file_name)
{
    //为在磁盘上保存connectome建立文件夹
    std::string pwd = "./mid_data/";
    if (!boost::filesystem::exists(pwd))
    {
        boost::filesystem::create_directories(pwd);
    }

    std::ofstream file(pwd + file_name + ".dat", std::ios_base::trunc);
    if (file.is_open() == false)
    {
        cout << "Error Matrix::dump : fail to open " + file_name << endl;
        exit(-1);
    }

    boost::archive::text_oarchive oa(file);
    oa << this->data;
    file.close();

    // cout<<"Dumping "<<file_name<<":"<<endl;

    // free connectome
    map<pair<int, int>, float, compare_col_first> temp = {};
    temp.swap(this->data);
}

void Matrix::load(string file_name)
{
    std::ifstream file("./mid_data/" + file_name + ".dat", std::ios_base::in);
    if (file.is_open() == false)
    {
        cout << "Error Matrix::load : fail to open " + file_name << endl;
        exit(-1);
    }

    boost::archive::text_iarchive ia(file);
    this->data.clear();
    ia >> this->data;

    file.close();
    // cout<<"Loading "<<file_name<<":"<<endl;
    // print_mat(connectome);
}

void Matrix::print()
{
    for (const auto &it : data)
    {
        cout << "((" << it.first.first << ", " << it.first.second << "), " << it.second << ") " << endl;
    }
}

// do not check index, if not exist, return 0
float Matrix::read_at(int row, int col)
{
    auto it = this->data.find(pair<int, int>(row, col));
    if (it != this->data.end())
    {
        return it->second;
    }
    else
    {
        return 0;
    }
}

void Matrix::clear()
{
    this->data.clear();
    this->max_row = 0;
    this->max_col = 0;
}

void Matrix::set_shape(int row, int col)
{
    this->max_row = row;
    this->max_col = col;
}

// write if there is element
void Matrix::multiply_at(int row, int col, const float &rate)
{
    auto it = this->data.find(pair<int, int>(row, col));
    if (it != this->data.end())
    {
        (*it).second *= rate;
    }
}

/**
 * @brief sum data[row][col], for row in winners
 */
float Matrix::get_col_winner_sum(int col, flat_set<int> &winners) //把winners 换成unordered_set?
{
    float sum = 0;
    auto upper_bound_it = this->data.upper_bound(pair<int, int>(-1, col + 1));
    for (auto it = this->data.lower_bound(pair<int, int>(0, col)); it != upper_bound_it; ++it)
    {
        if (winners.find(it->first.first) != winners.end()) // if it->first.first (row) is in winners
        {
            sum += it->second;
        }
    }
    return sum;
}

/*
void dump_connectome(std::string file_name, Matrix &connectome)
{
    std::ofstream file("./mid_data/" + file_name + ".dat", std::ios_base::trunc);
    if (file.is_open() == false)
    {
        cout << "Error dump_connectome: fail to open " + file_name << endl;
        exit(-1);
    }

    boost::archive::text_oarchive oa(file);
    oa << connectome;
    file.close();

    // cout<<"Dumping "<<file_name<<":"<<endl;
    // print_mat(connectome);

    // free connectome
    Matrix temp = {};
    temp.swap(connectome);
}

void load_connectome(std::string file_name, Matrix &connectome)
{
    std::ifstream file("./mid_data/" + file_name + ".dat", std::ios_base::in);
    if (file.is_open() == false)
    {
        cout << "Error load_connectome: fail to open " + file_name << endl;
        exit(-1);
    }

    boost::archive::text_iarchive ia(file);
    ia >> connectome;

    file.close();
    // cout<<"Loading "<<file_name<<":"<<endl;
    // print_mat(connectome);
}
*/

/**
 * @brief samples the truncated Normal PDF.
 *  Modified from: truncated_normal_ab_sample
 * @param size the required size of sample vector
 * @param start_index the start index of potential_new_winners
 * @note 可以考虑把double 改成float
 */
inline void truncated_normal_ab_sample_vec(deque<pair<int, float>> &potential_new_winners, double mu, double sigma, double a, double b,
                                           std::mt19937 &generator, int size, int start_index)
{
    double alpha;
    double alpha_cdf;
    double beta;
    double beta_cdf;
    double u;
    double x;
    double xi;
    double xi_cdf;
    int winner_size = 0;

    alpha = (a - mu) / sigma;
    beta = (b - mu) / sigma;

    alpha_cdf = normal_01_cdf(alpha);
    beta_cdf = normal_01_cdf(beta);

    for (int i = 0; i < size; ++i)
    {
        u = r8_uniform_01(generator());
        xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf);
        xi = normal_01_cdf_inv(xi_cdf);
        potential_new_winners.push_back(pair<int, float>(start_index + i, mu + sigma * xi));
    }
}

template <class T>
string print_vec(const T &vec)
{
    string output = "";
    for (const auto &it : vec)
    {
        output += std::to_string(it) + ",";
    }
    output += ";";
    return output;
}

template <class T>
void print_mat(const unordered_map<int, deque<T>> &mat)
{
    for (const auto &vec : mat)
    {
        // cout<<"Debug print_mat:vec.size="<<vec.size()<<endl;
        cout << print_vec(vec.second) + "\n";
    }
}

string print_vec_str(const deque<string> &vec)
{
    string output = "";
    for (const auto &it : vec)
    {
        output += "\"" + it + "\",";
    }
    // output += "\n";
    return output;
}

template <class T>
string print_vec_pair(const deque<pair<int, T>> &vec)
{
    string output = "";
    for (const auto &it : vec)
    {
        output += std::to_string(it.first) + ":" + std::to_string(it.second) + ",";
    }
    // output += "\n";
    return output;
}

inline bool compare_potential_winners(pair<int, float> p1, pair<int, float> p2)
{
    return p1.second > p2.second;
}

/**
 * @brief generate a vector of non-repeating random integers range in [0,range), range >= size must be established when calling
 * @param index   the vector to place the result
 * @param range     the range of ramdom numbers
 * @param size      the required size of the index
 */
inline void gen_random_index(std::vector<int> &index, int range, int size)
{
    if (size == 0 || range == 0)
    {
        // std::cerr << "Warning gen_random_index: range=" << range << " ,size=" << size << endl;
        return;
    }

    if (size < 0 || range - 1 < 0)
    {
        std::cerr << "Error gen_random_index: range=" << range << " ,size=" << size << endl;
        exit(-1);
    }

    index.clear();
    index.reserve(size);
    //蓄水池取样 可以使用Algorithm L https://en.wikipedia.org/wiki/Reservoir_sampling 再进行优化
    //另外，在之后以p为概率的采样中可以使用%，而不用distrubution： 不好用，因为p是浮点数而不是分数
    for (int i = 0; i < size; ++i)
    {
        index.push_back(i);
    }
    int j = 0;
    for (int i = size; i < range; ++i)
    {
        // Pick a random index from 0 to i.
        j = generator() % (i + 1);

        // If the randomly picked index is smaller than size,
        // then replace the element present at the index
        if (j < size)
        {
            index[j] = i;
        }
    }
}

/**
 * @brief generate a sample of vec_in with size elements
 * @param set_out    the set to place the result
 * @param set_in     the set to sample
 * @param size       size of the required sample
 */
inline void gen_random_sample(unordered_set<int> &set_out, const flat_set<int> &vec_in, long unsigned int size)
{
    vector<int> sample_index_set;
    if (size > vec_in.size())
    {
        std::cerr << "Error gen_random_sample: size must smaller than vector length. size=" << size << endl;
        return;
    }
    if (size == 0)
    {
        // std::cerr << "Warning gen_random_sample: size == 0" << endl;
        return;
    }
    if (size < 0)
    {
        std::cerr << "Error gen_random_sample: size=" << size << endl;
        exit(-1);
    }

    gen_random_index(sample_index_set, vec_in.size(), size);
    // cout<<"sample_index_set:"<<print_vec(sample_index_set)<<endl;

    set_out.clear();
    auto it = vec_in.begin();
    for (const auto &index : sample_index_set)
    {
        set_out.insert(*(it + index)); //这样写会不会越界？
        //对于winners，在这个随机取样函数中需要用下标取样（或其他随机取样方法），而在Matrix的get_col_winner_sum中需要快速查找
        //要不用vector排序后用二分查找？
        //或者用unordered_set挨个bucket取，取够为止？
        //用boost::container::flat_set
    }
}

/**
 * @brief padding by adding zeros at the end of each row
 * @param mat       the matrix need to be padded
 * @param new_size  the number of zeros need to be added for each row + the number of origin zeros of each row
 */
inline void padding_row(unordered_map<int, deque<float>> &mat, int new_size)
{
    if (new_size < mat.size())
    {
        std::cout << "Warning : padding_row new_size < mat.size():" << new_size << "<" << mat.size() << endl;
        //如果没有出错，应该只有在第一轮中会出现这样的现象
        return;
    }

    deque<float> null_vec(mat[0].size(), 0);
    // mat.reserve(mat.size() + pad_size);
    for (int i = mat.size(); i < new_size; ++i)
    {
        mat[i] = null_vec;
    }
}

/**
 * @brief padding by adding zeros for extra pad_size rows
 * @param mat       the matrix to be padded
 * @param new_size  the number of rows to be added + the number of original rows
 */
inline void padding_col(unordered_map<int, deque<float>> &mat, int new_size)
{
    if (new_size < mat[0].size())
    {
        std::cout << "Warning: padding_col new_size < mat[0].size():" << new_size << "<" << mat[0].size() << endl;
        // exit(-1);
        return;
    }
    // cout<<"before padding_col: mat[0].size="<<mat[0].size()<<endl;

    // cout<<"padding_col: pad_size="<<pad_size<<endl;
    // mat.reserve(mat.size() + pad_size);

    //不知道为什么，不能用(auto it: mat),否则可以用it正确写，但不能用mat访问
    //应该是auto默认是const，要加上&
    for (auto &it : mat)
    {
        // cout << "padding_col: padding mat[" << it->first << "]" << endl;
        it.second.resize(new_size, 0); // 前面已经保证pad_size >= it->second.size()
    }

    // cout<<"after padding_col: mat[0].size="<<mat[0].size()<<endl;
}

Stimulus::Stimulus(int k0)
{
    this->k = k0;
}

Area::Area(string name0, int n0, int k0, float beta0)
{
    this->name = name0; // deep copy with (copy on write) optimization
    this->n = n0;
    this->k = k0;
    this->beta = beta0;

    this->w = 0;
    this->saved_w = {0};

    this->new_w = 0;

    this->num_first_winners = -1;
    this->fixed_assembly = false;
    this->explicit_ = false;
}
void Area::update_winners()
{
    this->winners.swap(this->new_winners);
    if (this->explicit_ == false)
    {
        this->w = this->new_w;
    }
}
void Area::update_stimulus_beta(string name0, float new_beta0)
{
    this->stimulus_beta.insert(pair<string, float>(name0, new_beta0)); // insert 会自动检查键是否已经存在，如果存在则插入失败
}
void Area::update_area_beta(string name0, float new_beta0)
{
    this->area_beta.insert(pair<string, float>(name0, new_beta0));
}
void Area::fix_assembly()
{
    if (this->winners.size() == 0)
    {
        std::cerr << "Error: Area " + this->name + " does not have assembly, cannot fix." << endl;
        return;
    }
    this->fixed_assembly = true;
}
void Area::unfix_assembly()
{
    this->fixed_assembly = false;
}

Brain::Brain(float p0, bool save_size0, bool save_winners0)
{
    this->p = p0;
    this->save_size = save_size0;
    this->save_winners = save_winners0;
    this->no_plasticity = false;
}

/**
 * @brief generate a deque of random integers in the range [0,k], where each value represents the number of successes in a sequence of k trials (each with a probability of success equal to p).
 * @param vec_out the deque to put results
 * @param k  the number of trials
 * @param p  the probability of success of each trial
 * @param size the size of the deque
 */
void gen_binomial_vec(deque<float> &vec_out, int k, float p, int size)
{
    static std::binomial_distribution<int> binomial_distrib(k, p); //为什么不用正态分布？

    vec_out.clear();
    // vec_out.reserve(size);
    for (int i = 0; i < size; ++i)
    {
        vec_out.push_back(float(binomial_distrib(generator))); // emplace_back() vs push_back() https://developer.aliyun.com/article/771502
    }
}

/**
 * @brief generate a matrix of random integers in the range [0,k], where each value represents the number of successes in a sequence of k trials (each with a probability of success equal to p).
 * @param mat_out the matrix to put results
 * @param k  the number of trials
 * @param p  the probability of success of each trial
 * @param size the size of the matrix in (row,col)
 */
void gen_binomial_mat(Matrix &mat_out, int k, float p, pair<int, int> size)
{
    static std::binomial_distribution<int> binomial_distrib(k, p); //为什么不用正态分布？

    deque<float> null_vec;
    int row = size.first;
    int col = size.second;

    mat_out.clear();
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            mat_out.at(i, j) = binomial_distrib(generator);
        }
    }
}

// # Add stimulus to brain. The weights of the connections to an explicit area of the brain are set randomly
// # according to a binomial distribution with parameters k and the value p of the brain. The plasticity of
// # these connections are set equal to the default plasticity of the area.这里说的应该是beta
// did not examine whether the name0 is used （之后加上）
void Brain::add_stimulus(string name0, int k0) // k0是刺激区域在本轮中，激活的神经元数量
{
    this->stimuli[name0] = Stimulus(k0); // dequey应该不用我来维护Stimulus
    unordered_map<string, deque<float>> new_connectomes;
    deque<float> null_vec;
    string key;
    Area *area_ptr;
    for (auto &areas_it : this->areas)
    {
        key = areas_it.first;
        area_ptr = &(areas_it.second);

        if (area_ptr->explicit_ == true)
        {
            gen_binomial_vec(new_connectomes[key], k0, this->p, area_ptr->n);
            // cout << "Debug: new_connectomes[" + key + "]=" + print_vec(new_connectomes[key]) << endl;
        }
        else
        {
            // null_vec.resize(area_ptr->k, 1); //这里k的位置其实是填各area已有记录的神经元数量，不妨先填为k
            new_connectomes[key] = null_vec;
        }
        // end if explicit

        area_ptr->stimulus_beta[name0] = area_ptr->beta; //用下标的方式修改unordered_map，如果已有该键会直接覆盖。如果用insert且已有改键则会报错。
    }                                                    // end for areas

    this->stimuli_connectomes[name0] = new_connectomes; // stimuli_connectomes[stim][area]
}

void generate_null_mat(unordered_map<int, deque<float>> &null_mat, int num_of_row, int num_of_col)
{
    null_mat.clear();
    deque<float> null_vec(num_of_col, 0); //初始化为0
    for (int i = 0; i < num_of_row; ++i)
    {
        null_mat[i] = null_vec;
    }
}

// # Add a non-explicit area to the brain. Since the area is not explicit, all connections from the
// # stimuli of the brain and from/to all the areas of the brain are initially set to empty (they
// # will be set during the project phase in order to improve performance).
void Brain::add_area(string name0, int n0, int k0, float beta0)
{
    this->areas[name0] = Area(name0, n0, k0, beta0);

    deque<float> null_vec;

    for (auto &stimuli_connectomes_it : this->stimuli_connectomes)
    {
        stimuli_connectomes_it.second[name0] = null_vec; //初始化为空
        this->areas[name0].stimulus_beta[stimuli_connectomes_it.first] = beta0;
    }

    unordered_map<string, Matrix> new_connectomes; // connectomes[name0][other_area]
    int other_area_size = 0;
    Matrix null_mat; //会出现空指针，需要初始化
    string key;
    Area *area_ptr;
    for (auto &areas_it : this->areas)
    {
        key = areas_it.first;
        area_ptr = &(areas_it.second);

        other_area_size = area_ptr->k; //原版是other_area_size = 0
        if (area_ptr->explicit_ == true)
        {
            other_area_size = area_ptr->n;
        }

        null_mat.set_shape(0, other_area_size); // np.zeros((0, other_area_size)) 用于给后面的padding提供row或col的参数
        // null_mat.dump("connectome[" + name0 + "][" + key + "]");
        new_connectomes[key] = null_mat; //包括了自环的情况

        if (key != name0)
        {
            null_mat.set_shape(other_area_size, 0); // np.zeros((other_area_size, 0))
            // null_mat.dump("connectome[" + key + "][" + name0 + "]");
            this->connectomes[key][name0] = null_mat; //(other_area_size, 0) 会在之后的每次新增中用padding补充并初始化成(other_area_size, w)
        }

        area_ptr->area_beta[name0] = area_ptr->beta;
        this->areas[name0].area_beta[key] = beta0;
    } // end for areas

    this->connectomes[name0] = new_connectomes;
}

// # Add an explicit area to the brain. Since the area is explicit, the weights of all connections
// # from a stimuli of the brain the new area are initially set randomly according to a
// # binomial distribution with parameters the value k of the stimulus and the value p of the brain
// # (with the default plasticity). The weights of all connections from/to an explicit area of the brain
// # to the new area are initially and fully set randomly according to a binomial distribution
// # with parameters 1 and the value p of the brain. The weights of all connections from/to a
// # non-explicit area of the brain to the new area are initially set to empty.
// # In all cases, the plasticity of the connections is set to the default plasticity.
// # The number of winners of the new area is set equal to the number of its neurons.
void Brain::add_explicit_area(string name0, int n0, int k0, int beta0)
{
    this->areas[name0] = Area(name0, n0, k0, beta0);
    this->areas[name0].explicit_ = true;

    string stim_name;
    deque<float> *connectome_ptr;
    for (auto &stimuli_connectomes_it : this->stimuli_connectomes)
    {
        for (auto &stimuli_it : stimuli_connectomes_it.second)
        {
            stim_name = stimuli_it.first;
            connectome_ptr = &(stimuli_it.second);

            gen_binomial_vec(*connectome_ptr, this->stimuli[stim_name].k, this->p, n0); //这里的stimuli[].k只是随机生成用的参数，不是规模
            this->areas[name0].stimulus_beta[stim_name] = beta0;
        }
    }

    string key;
    Area *area_ptr;
    int other_n;
    Matrix null_mat;
    unordered_map<string, Matrix> new_connectomes;
    for (auto &areas_it : this->areas)
    {
        key = areas_it.first;
        area_ptr = &(areas_it.second);
        if (key == name0)
        {
            gen_binomial_mat(new_connectomes[key], 1, p, pair<int, int>(n0, n0));
            // new_connectomes[key].dump("connectome[" + name0 + "][" + key + "]");
        }
        if (key != name0)
        {
            if (area_ptr->explicit_ == true)
            {
                other_n = area_ptr->n;
                gen_binomial_mat(new_connectomes[key], 1, p, pair<int, int>(n0, other_n));
                // new_connectomes[key].dump("connectome[" + name0 + "][" + key + "]");

                gen_binomial_mat(this->connectomes[key][name0], 1, p, pair<int, int>(other_n, n0));
                // this->connectomes[key][name0].dump("connectome[" + key + "][" + name0 + "]");
            }
            else
            {
                new_connectomes[key] = null_mat;
                // new_connectomes[key].dump("connectome[" + name0 + "][" + key + "]");
            }

            this->connectomes[key][name0] = null_mat;
            // this->connectomes[key][name0].dump("connectome[" + key + "][" + name0 + "]");
        } // end if key != name0
        area_ptr->area_beta[name0] = area_ptr->beta;
        this->areas[name0].area_beta[key] = beta0;
    } // end for areas

    this->connectomes[name0] = new_connectomes;
    // Explicitly set w to n so that all computations involving this area are explicit.
    this->areas[name0].w = n0;
}

// # Update the plasticities of the connections between stimuli and areas. Each area update consists of
// # of the destination area and a list of update rules: each update rule specifies the source area
// # and the new plasticity. Each stimulus update consists of the destination area and a list of update rules:
// # each update rule specifies the source stimulus and the new plasticity.
//#area_update_unordered_map consists of area1 : list[(area2, new_beta)]
void Brain::update_plasticities(unordered_map<string, unordered_map<string, float>> area_update_unordered_map, unordered_map<string, unordered_map<string, float>> stim_update_unordered_map)
{
    // area_update_unordered_map consists of area1: list[ (area2, new_beta) ]
    // represents new plasticity FROM area2 INTO area1
    for (const auto &area_update_it : area_update_unordered_map)
    {
        for (const auto &rule_it : area_update_it.second)
        {
            this->areas[area_update_it.first].area_beta[rule_it.first] = rule_it.second;
        }
    }

    // stim_update_unordered_map consists of area: list[ (stim, new_beta) ]f
    // represents new plasticity FROM stim INTO area
    for (const auto &stim_update_it : stim_update_unordered_map)
    {
        for (const auto &rule_it : stim_update_it.second)
        {
            this->areas[stim_update_it.first].stimulus_beta[rule_it.first] = rule_it.second;
        }
    }
}

// # Execute the project from stimuli and/or areas to areas. For each stimulus (key) in the first dictionary,
// # the list (value) of areas to which the stimulus has the project is specified. For each area (key),
// # in the second dictionary, the list (value) of areas to which the area has the project is specified.
// # The function collects, for each area, the list of stimuli and areas that project to it (basically, it
// # computes the inverse of the input mappings). Then, for each area which has to "receive" a projection
// # (from either stimuli or areas), it invokes the function which actually performs the projection (this
// # function returns the number of winners in the destination area). If the new winners have to be later
// # analysed, then their list is appended to the the list of lists of winners of the area. When everything
// # has been done, the function updates the destination areas.
void Brain::project(const unordered_map<string, deque<string>> &stim_to_area, const unordered_map<string, deque<string>> &area_to_area, bool verbose)
{
    static std::binomial_distribution<int> binomial_distrib(1, this->p);

    //  Validate stim_area, area_area well defined
    //   stim_to_area: {"stim1":["A"], "stim2":["C","A"]}
    // 	 area_to_area: {"A":["A","B"],"C":["C","A"]}
    unordered_map<string, deque<string>> stim_in; //默认没有定义的键返回空的deque
    unordered_map<string, deque<string>> area_in; //默认没有定义的键返回空的deque

    string stim;
    string to_area;
    for (const auto &projection : stim_to_area) // auto 是默认const,需要加上&，除非是vector
    {
        stim = projection.first;

        if (stimuli.find(stim) == stimuli.end())
        { // if stim not in stimuli
            std::cerr << "Error: " + stim + " not in Brain.stimuli (stim_to_area)" << endl;
            exit(-1);
        }

        for (const auto &to_area : projection.second)
        {
            if (this->areas.find(to_area) == this->areas.end())
            { // if to_area not in areas
                std::cerr << "Error: " + to_area + " not in Brain.areas (stim_to_area)" << endl;
                exit(-1);
            }
            stim_in[to_area].push_back(stim); // append
        }
    } // end for stim_to_area

    string from_area;
    for (const auto &projection : area_to_area)
    {
        from_area = projection.first;

        if (this->areas.find(from_area) == this->areas.end())
        {
            std::cerr << "Error: " + from_area + " not in Brain.areas (area_to_area)" << endl;
            exit(-1);
        }

        for (const auto &to_area : projection.second)
        {
            if (areas.find(to_area) == areas.end())
            {
                std::cerr << "Error: " + to_area + " not in Brain.areas (area_to_area 2)" << endl;
                exit(-1);
            }
            area_in[to_area].push_back(from_area);
        }
    } // end for area_to_area

    set<string> to_update; //所有需要更新的area
    for (const auto &it : stim_in)
    {
        to_update.insert(it.first);
    }
    for (const auto &it : area_in)
    {
        to_update.insert(it.first);
    }

    int num_first_winners = 0;
    deque<int> *saved_winners_ptr;
    flat_set<int> *new_winners_ptr;
    // unordered_map<string, Matrix> temp_connectomes; // defalut as connectomes[other_area][update_area], omit [update_area]
    Matrix *temp_mat_ptr;
    Matrix null_mat;
    Area *from_area_ptr;
    Area *area0_ptr;
    int columns = 0;
    for (const auto &update_area : to_update)
    {
        // temp_connectomes.clear();

        for (const auto &from_area : area_in[update_area])
        {
            // 读入connectomes[from_area][name0]
            // temp_connectomes[from_area] = &(this->connectomes[from_area][to_update]);

            // 根据other_area.saved_w和other_area.w,写入上一轮的connectomes[name0][other_area] (padding_row)
            temp_mat_ptr = &(this->connectomes[from_area][update_area]);
            from_area_ptr = &(this->areas[from_area]);
            area0_ptr = &(this->areas[update_area]);

            // if (std::find(from_areas.begin(), from_areas.end(), from_area) == from_areas.end())
            if (temp_mat_ptr->get_max_col() < area0_ptr->w) //在第一次读取时会下标越界，需要更换Matrix实现
            {                                               // if other_area not in from_areas of update_area
                // padding_col(*temp_mat_ptr, area0_ptr->w);

                for (int j = 0; j < temp_mat_ptr->get_max_row(); ++j)
                {
                    for (int i = temp_mat_ptr->get_max_col(); i < area0_ptr->w; ++i)
                    {
                        // cout << "connectomes[" << other_area << "][" << name0 << "].at(" << j << ").at(" << i << ")" << endl;

                        // add num_first_winners rows, all bernoulli with probability p
                        // (*temp_mat_ptr).at(j).at(i) = binomial_distrib(generator); // np.random.binomial(1, self.p)
                        if (int(binomial_distrib(generator)) == 1)
                        {
                            temp_mat_ptr->at(j, i) = 1; // np.random.binomial(1, self.p)
                        }
                    }
                }

                temp_mat_ptr->set_col(area0_ptr->w);
            }

            if (temp_mat_ptr->get_max_row() < from_area_ptr->w)
            {
                // padding_row(*temp_mat_ptr, from_area_ptr->w);

                for (int i = temp_mat_ptr->get_max_row(); i < from_area_ptr->w; ++i)
                {
                    for (int j = 0; j < temp_mat_ptr->get_max_col(); ++j)
                    {
                        if (int(binomial_distrib(generator)) == 1)
                        {
                            temp_mat_ptr->at(i, j) = 1; // np.random.binomial(1, self.p);
                        }
                    }
                }

                temp_mat_ptr->set_row(from_area_ptr->w);
            }

            if (verbose == true)
            {
                cout << "Connectome of " + from_area + " to " + update_area + " is now:";
                // temp_connectomes[from_area].print();
            }
        }

        num_first_winners = this->project_into(this->areas[update_area], stim_in[update_area], area_in[update_area], verbose);
        // pthread 有读写锁，但如果把project_into分段执行可能就不用锁了
        //有很多C标准库函数不是线程安全的
        // std::thread 可以批量创建线程，加锁并且有thread_local

        /*
        //将对temp_connectomes的修改写回硬盘
        for (const auto &from_area : area_in[update_area])
        {
            // dump_connectome("connectome[" + from_area + "][" + update_area + "]", temp_connectomes[from_area]);
            temp_connectomes[from_area].dump("connectome[" + from_area + "][" + update_area + "]");
        }
        */

        this->areas[update_area].num_first_winners = num_first_winners;
        cout << "num_first_winners=" << num_first_winners << endl;

        if (this->save_winners == true)
        {
            saved_winners_ptr = &(this->areas[update_area].saved_winners);
            new_winners_ptr = &(this->areas[update_area].new_winners);

            saved_winners_ptr->insert(saved_winners_ptr->end(), new_winners_ptr->begin(), new_winners_ptr->end()); // append new_winners to saved_winners
        }
    } // end for to_update
    // once done everything, for each area in to_update: area.update_winners()
    for (const auto &update_area : to_update)
    {
        this->areas[update_area].update_winners();
        if (save_size == true)
        {
            this->areas[update_area].saved_w.push_back(this->areas[update_area].w);
        }
    } // end for to_update
}

int Brain::project_into(Area &area0, deque<string> from_stimuli, deque<string> from_areas, bool verbose)
{
    // projecting everything in from stim_in[area] and area_in[area]
    // calculate: inputs to self.connectomes[area] (previous winners)
    // calculate: potential new winners, Binomial(sum of in sizes, k-top)
    // k top of previous winners and potential new winners
    // if new winners > 0, redo connectome and intra_connectomes
    // have to wait to replace new_winners
    std::clog << "Projecting " + print_vec_str(from_stimuli) + " and " + print_vec_str(from_areas) + " into " + area0.name << endl;

    // If projecting from area with no assembly, complain.
    for (const auto &from_area : from_areas)
    {
        if ((this->areas[from_area].winners.size() == 0) || (this->areas[from_area].w == 0))
        {
            cout << "Projecting from area with no assembly: " + from_area << endl;
        }
    }

    string name0 = area0.name;
    deque<pair<int, float>> prev_winner_inputs(area0.w); // size=area0.w init to 0.0,可以尝试使用unordered_map<int,float>
    for (int i = 0; i < area0.w; ++i)
    {
        prev_winner_inputs[i].first = i;
        prev_winner_inputs[i].second = 0;
    }

    deque<float> *stim_inputs_ptr;
    for (const auto &stim : from_stimuli)
    {
        stim_inputs_ptr = &(this->stimuli_connectomes[stim][name0]); //这里需要读取stimuli
        // cout << "Debug: area0.w=" << area0.w << " ,stim_inputs.size=" << stim_inputs_ptr->size() << " ,stimuli_connectomes[" + stim + "][" + name0 + "]" << endl;
        for (int i = 0; i < area0.w; ++i)
        {
            prev_winner_inputs[i].second += (*stim_inputs_ptr).at(i);
            // if ((*stim_inputs_ptr).at(i) != 0)
            //     cout << "prev_winner_inputs[" << i << "].second +=" << (*stim_inputs_ptr).at(i) << " stimulus at(" << i << ")" << endl;
        }
    } // end for from_stimuli
    Matrix *connectome0_ptr;
    for (const auto &from_area : from_areas)
    {
        connectome0_ptr = &(this->connectomes[from_area][name0]); // two-dim float
        // cout << "connectomes[" << from_area << "][" << name0 << "]=(" << connectomes[from_area][name0].size() << "," << connectomes[from_area][name0][0].size() << ")" << endl;

        for (int i = 0; i < area0.w; ++i)
        {
            prev_winner_inputs[i].second += connectome0_ptr->get_col_winner_sum(i, this->areas[from_area].winners);
            // if (connectome0_ptr->get_col_winner_sum(i, this->areas[from_area].winners) != 0)
            //         cout << "prev_winner_inputs[" << i << "].second +=" << connectome0_ptr->get_col_winner_sum(i, this->areas[from_area].winners) << " connectome at col "<< i<< endl;
        }

    } // end for from_areas
    if (verbose == true)
    {
        cout << "prev_winner_inputs: " << print_vec_pair(prev_winner_inputs) << endl;
    }

    // simulate area.k potential new winners if the area is not explicit
    int total_k = 0;
    deque<int> input_sizes;
    int num_inputs = 0;
    int effective_k = 0;
    int effective_n = 0;
    float alpha = 0;

    float std = 0;
    float mu = 0;
    // float a = 0;
    // float b = 0;
    deque<pair<int, float>> potential_new_winners; // <index,input weight>,size=area0.k,可能会在内存放不下 , directly push into prev_winner_inputs

    deque<pair<int, float>> *all_potential_winners_ptr = NULL;

    if (area0.explicit_ == false)
    {
        // # Compute the number of input stimuli and areas, the total number of input connectomes,
        // # and the number of input connectomes for each input stimulus and area.

        for (const auto &stim : from_stimuli)
        {
            total_k += this->stimuli[stim].k; //这里应该是记录所有刺激中产生了信号的神经元数量
            input_sizes.push_back(this->stimuli[stim].k);
            num_inputs += 1;
        } // end for from_stimuli

        for (const auto &from_area : from_areas)
        {
            effective_k = this->areas[from_area].winners.size();
            total_k += effective_k; //记录所有脑区范围内本轮被激活的神经元的数量
            input_sizes.push_back(effective_k);
            num_inputs += 1;
        } // end for from_areas

        if (verbose == true)
        {
            // total_k 是全脑本轮被激活的神经元总数，input_sizes是分别统计所有刺激和所有脑区本轮活动的神经元
            cout << "total_k = " + std::to_string(total_k) + " and input_sizes = " + print_vec(input_sizes) << endl;
        }

        // # Compute the potential new k winners of the area to which we are going to project.
        // # To this aim compute the threshold alpha for inputs that are above (n-k)/n percentile,
        // # use the normal approximation, between alpha and total_k, round to integer, and
        // # create the k potential_new_winners.
        // function compute_potential_new_winners(b::Brain, to_area::Area, total_k::Int)

        effective_n = area0.n - area0.w; //还没有被记录的神经元数量
        // Threshold for inputs that are above (n-k)/n percentile.
        // self.p can be changed to have a custom connectivity into this brain area.
        // alpha = binom.ppf((float(effective_n-area0.k)/effective_n), total_k, self.p)
        // The.ppf () function calculates the probability for a given normal distribution value
        //, while the.cdf () function calculates the normal distribution value for which a given probability is the required value. These are inverse of each other in this particular sense.
        alpha = quantile(binom_round_up(total_k, this->p), float(effective_n - area0.k) / effective_n);
        //这样求的alpha有问题
        if (verbose == true)
        {
            cout << "Alpha = " << alpha << endl;
        }

        // use normal approximation, between alpha and total_k, round to integer
        // create k potential_new_winners
        std = std::sqrt(total_k * (this->p) * (1.0 - (this->p)));
        mu = total_k * (this->p);
        // a = float(alpha - mu) / std;
        // b = float(total_k - mu) / std;
        static std::normal_distribution<float> normal_distrib(mu, std);

        // 以下在实现python中的potential_new_winners = truncnorm.rvs(a, b, scale=std, size=area.k)
        //多次对标准正态分布采样，保留在区间[a,b]中的采样点,返回浮点数
        // cout<<"mu="<<mu<<" ,std="<<std<<" ,alpha="<<alpha<<" ,total_k="<<total_k<<endl;
        // for (int i = 0; i < area0.k; ++i)
        // {
        //     temp_winner = float(truncated_normal_ab_sample(mu, std, alpha, total_k, generator())); //这里还可以改一下源码，让它直接生成一个数组
        //     // cout<<temp_winner<<" ";
        //     potential_new_winners.push_back(pair<int, float>(area0.w + winner_size, temp_winner));
        // }
        truncated_normal_ab_sample_vec(potential_new_winners, mu, std, alpha, total_k, generator, area0.k, area0.w);

        if (verbose == true)
        {
            cout << "potential_new_winners: " << print_vec_pair(potential_new_winners) << endl;
        }

        // take max among prev_winner_inputs, potential_new_winners
        // get num_first_winners (think something small)
        // can generate area.new_winners, note the new indices
        // all_potential_winners = prev_winner_inputs + potential_new_winners
        // pair<int,float>
        // cout << "prev_winner_inputs" << print_vec_pair(prev_winner_inputs) << "\npotential_new_winners:" << print_vec_pair(potential_new_winners) << endl;

        // prev_winner_inputs.reserve(prev_winner_inputs.size() + potential_new_winners.size()); //预先分配空间，加快复制速度
        // for (const auto &it : potential_new_winners)
        // {
        //     prev_winner_inputs.push_back(it);
        // }
        prev_winner_inputs.insert(prev_winner_inputs.end(), potential_new_winners.begin(), potential_new_winners.end()); //这一步应该有更好的合并方式

    } // end if area0 not explicit
    all_potential_winners_ptr = &prev_winner_inputs;

    // 生成一组下标，把下标按照all_potential_winners中的值进行排序,只排出top area0.k
    // new_winner_indices = heapq.nlargest(area.k, range(len(all_potential_winners)), all_potential_winners.__getitem__) #为什么要的结果是下标？all_potential_winners里的元素顺序是对应connectome里的下标吗？
    std::partial_sort(all_potential_winners_ptr->begin(), all_potential_winners_ptr->begin() + area0.k, all_potential_winners_ptr->end(), compare_potential_winners);
    // cout<<print_vec_pair(*all_potential_winners_ptr);

    vector<int> new_winner_indices;
    new_winner_indices.reserve(area0.k); //预先分配空间
    for (int i = 0; i < area0.k; ++i)
    {
        new_winner_indices.push_back((*all_potential_winners_ptr)[i].first);
    }
    // cout << "new_winner_indices.size=" << new_winner_indices.size() << endl;
    // cout<<"new_winner_indices:\n"<<print_vec(new_winner_indices)<<endl;

    int num_first_winners = 0;
    deque<int> first_winner_inputs;
    if (area0.explicit_ == false)
    {
        for (auto &new_winner_index : new_winner_indices)
        {
            if (new_winner_index >= area0.w)
            {
                first_winner_inputs.push_back(potential_new_winners[new_winner_index - area0.w].second);
                new_winner_index = area0.w + num_first_winners;
                num_first_winners++;
            }
        }
        /*
        for (int i = 0; i < area0.k; ++i)
        {
            if (new_winner_indices[i] >= area0.w)
            {
                first_winner_inputs.push_back(potential_new_winners[new_winner_indices[i] - area0.w].second);
                new_winner_indices[i] = area0.w + num_first_winners;
                num_first_winners++;
            }
        }
        */
    } // end if not explicit

    // cout << "len(first_winner_inputs)=" << first_winner_inputs.size() << endl;
    // area0.new_winners = new_winner_indices; //如果用swap，由于new_winner_indices是函数内的变量，会出错
    area0.new_winners.clear();
    area0.new_winners.reserve(area0.k);
    area0.new_winners.insert(new_winner_indices.begin(), new_winner_indices.end());

    area0.new_w = area0.w + num_first_winners;

    // For experiments with a "fixed" assembly in some area.
    if (area0.fixed_assembly == true)
    {
        area0.new_winners.swap(area0.winners);
        area0.new_w = area0.w;
        first_winner_inputs.clear();
        num_first_winners = 0;
    }

    if (verbose == true)
    {
        cout << "new_winners: " << print_vec(area0.new_winners) << endl;
    }
    // for i in num_first_winners
    // generate where input came from
    // 	1) can sample input from array of size total_k, use ranges
    // 	2) can use stars/stripes method: if m total inputs, sample (m-1) out of total_k
    map<int, deque<int>> first_winner_to_inputs;
    // if (num_first_winners > 0)
    // {
    //     first_winner_to_inputs.resize(num_first_winners);
    // }

    vector<int> input_indices;
    int total_so_far = 0;
    int sum_proper_w = 0; // the number of w meeting certain criteria
    deque<int> inputs(num_inputs, 0);
    for (int i = 0; i < num_first_winners; ++i)
    {
        gen_random_index(input_indices, total_k, first_winner_inputs[i]); // 从range[0,total_k)中不重复随机取样first_winner_inputs[i]个
        inputs.clear();
        inputs.resize(num_inputs, 0);
        total_so_far = 0;
        // cout << "num_inputs=" << num_inputs << " num_first_winners=" << num_first_winners << endl;
        for (int j = 0; j < num_inputs; ++j)
        {
            sum_proper_w = 0;
            for (const auto &w : input_indices)
            {
                if ((total_so_far + input_sizes[j]) > w && w >= total_so_far)
                {
                    sum_proper_w++;
                }
            }
            inputs[j] = sum_proper_w;
            total_so_far += input_sizes[j];
        } // end for j in range(num_inputs)
        first_winner_to_inputs[i] = inputs;

        if (verbose == true)
        {
            cout << "For first_winner # " << i << " with input " << first_winner_inputs[i] << " split as so: \n"
                 << print_vec(inputs) << endl;
        }
    } // end for i < num_first_winners

    // connectome for each stim->area
    // add num_first_winners cells, sampled input * (1+beta)
    // for i in repeat_winners, stimulus_inputs[i] *= (1+beta)
    int m = 0;
    float stim_to_area_beta = 0;
    for (const auto &stim : from_stimuli)
    {
        if (num_first_winners > 0)
        {
            this->stimuli_connectomes[stim][name0].resize(area0.w + num_first_winners);
        }
        for (int i = 0; i < num_first_winners; ++i)
        {
            this->stimuli_connectomes[stim][name0].at(area0.w + i) = first_winner_to_inputs[i][m];
        }

        stim_to_area_beta = area0.stimulus_beta[stim];

        if (this->no_plasticity == true)
        {
            stim_to_area_beta = 0;
        }

        for (const auto &i : area0.new_winners)
        {
            this->stimuli_connectomes[stim][name0].at(i) *= (1 + stim_to_area_beta);
        }

        if (verbose == true)
        {
            cout << "stimuli_connectomes[" << stim << "][" << name0 << "]=" << print_vec(this->stimuli_connectomes[stim][name0]) << endl;
        }

        m++; //记录所有输入的数量，包括stimuli和area
    }        // end for stim in from_stimuli

    int from_area_w = 0;
    flat_set<int> *from_area_winners_ptr;
    int total_in = 0;
    unordered_set<int> sample_indices; //这里应该使用set或vector，用deque会在下面判断元素是否存在时出错
    static std::binomial_distribution<int> binomial_distrib(1, this->p);
    float area_to_area_beta_rate = 1;
    for (const auto &from_area : from_areas)
    {
        from_area_w = this->areas[from_area].w;

        from_area_winners_ptr = &(this->areas[from_area].winners); //这里需要读取其他area的winner，只读不写，造成数据冲突

        // padding col
        for (int i = 0; i < num_first_winners; ++i)
        {
            // cout << "first_winner_to_inputs[" << i << "][" << m << "]=" << first_winner_to_inputs[i][m] << endl;
            total_in = first_winner_to_inputs[i][m];

            // generate a sample of from_area_winners with number of total_in elements
            gen_random_sample(sample_indices, this->areas[from_area].winners, total_in);
            // cout << "sample_indices:" << print_vec(sample_indices) << endl;

            // cout << "connectomes[" << from_area << "][" << name0 << "]=(" << connectomes[from_area][name0].size() << "," << connectomes[from_area][name0][0].size() << ")" << endl;
            for (int j = 0; j < this->connectomes[from_area][name0].get_max_row(); ++j) // from_area_w
            {
                if (from_area_winners_ptr->find(j) == from_area_winners_ptr->end())
                { // if j not in from_area_winners
                    // cout << "connectomes[" << from_area << "][" << name0 << "].at(" << j << ").at(" << area0.w + i << ")" << endl;
                    if (binomial_distrib(generator) == 1)
                    {
                        this->connectomes[from_area][name0].at(j, area0.w + i) = 1.0; // np.random.binomial(1, self.p)
                    }
                }
                else if (sample_indices.find(j) != sample_indices.end())
                { // if j in sample_indices
                    // cout << "connectomes[" << from_area << "][" << name0 << "].at(" << j << ").at(" << area0.w + i << ")" << endl;
                    this->connectomes[from_area][name0].at(j, area0.w + i) = 1.0;
                }

            } // end for j
        }     // end for i
        this->connectomes[from_area][name0].set_col(area0.w + num_first_winners);

        area_to_area_beta_rate = area0.area_beta[from_area] + 1;
        if (this->no_plasticity == true)
        {
            area_to_area_beta_rate = 1.0;
        }

        for (const auto &i : area0.new_winners)
        {
            for (const auto &j : (*from_area_winners_ptr))
            {
                // cout << "connectomes[" << from_area << "][" << name0 << "].at(" << j << ").at(" << i << ")" << endl;
                this->connectomes[from_area][name0].multiply_at(j, i, area_to_area_beta_rate); // multiply for non 0 elements
            }
        }
        if (verbose == true)
        {
            // cout << "Connectome of " + from_area + " to " + name0 + " is now:\n";
            this->connectomes[from_area][name0].print();
        }

        m++;
    } // end for from_area in from_areas

    return num_first_winners;
}

#endif