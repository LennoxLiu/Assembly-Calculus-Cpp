#include "brain.hpp"
#include <iostream>
#include <ctime>

void project_sim(Brain &b, int n0 = 1000000, int k0 = 1000, float p0 = 0.01, float beta0 = 0.05, int t0 = 50, bool verbose = false)
{
    cout << "n=" << n0 << ",k=" << k0 << ",p=" << p0 << ",beta=" << beta0 << ",t=" << t0 << endl;

    b.add_stimulus("stim", k0);
    b.add_area("A", n0, k0, beta0);

    cout << "turn 0" << endl;
    b.project({
                  {"stim", {"A"}},
              },
              {}, verbose);

    for (int i = 0; i < t0 - 1; ++i)
    {
        cout << "turn" << i + 1 << endl;
        b.project({
                      {"stim", {"A"}},
                  },
                  {
                      {"A", {"A"}},
                  },
                  verbose);
    }
    // return b.areas["A"].saved_w;
}

void merge_sim(Brain &b, int n0 = 100000, int k0 = 317, float p0 = 0.01, float beta0 = 0.05, int max_t0 = 50, bool verbose = false)
{
    cout << "n=" << n0 << ",k=" << k0 << ",p=" << p0 << ",beta=" << beta0 << ",max_t=" << max_t0 << endl;
    b.add_stimulus("stimA", k0);
    b.add_stimulus("stimB", k0);
    b.add_area("A", n0, k0, beta0);
    b.add_area("B", n0, k0, beta0);
    b.add_area("C", n0, k0, beta0);

    unordered_map<string, deque<string>> stim_to_area;
    unordered_map<string, deque<string>> area_to_area;
    b.project({
                  {"stimA", {"A"}},
              },
              {},verbose);
    b.project({
                  {"stimB", {"B"}},
              },
              {},verbose);
    b.project({
                  {"stimA", {"A"}},
                  {"stimB", {"B"}},
              },
              {{"A", {"A", "C"}}, {"B", {"B", "C"}}},verbose);
    b.project({
                  {"stimA", {"A"}},
                  {"stimB", {"B"}},
              },
              {
                  {"A", {"A", "C"}},
                  {"B", {"B", "C"}},
                  {"C", {"C", "A", "B"}},
              },verbose);
    for (int i = 0; i < max_t0 - 1; ++i)
    {
        cout << "turn " << i << endl;
        b.project({
                      {"stimA", {"A"}},
                      {"stimB", {"B"}},
                  },
                  {
                      {"A", {"A", "C"}},
                      {"B", {"B", "C"}},
                      {"C", {"C", "A", "B"}},
                  },verbose);
    }
}

int main(int argc, char *argv[])
{
    clock_t start_time = clock();

    Brain *b = new Brain(0.01);
    // project_sim(*b, 100, 10, 0.01, 0.05, 10, false);
    merge_sim(*b,100000,317,0.01,0.05,50,false);

    cout << "Time overall:" << (clock() - start_time) / (double)CLOCKS_PER_SEC << "s" << endl;
    delete (b);
}