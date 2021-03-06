#ifndef PATHFINDER_H
#define PATHFINDER_H

#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <cstdlib>
#include <limits>

#include "main.hpp"
#include "transition_matrix.hpp"
#include "progressbar.hpp"

//! Paths represented by arrays and their probabilty as float
typedef std::map<std::vector<State>, double> Paths;
typedef std::map<std::vector<State>, double> Waiting_times;
typedef std::multimap<double, std::vector<State>> Paths_sorted;
class Pathfinder;  // forward declaration
typedef std::map<State, Pathfinder*> Pathfinder_map;

void simplify_path(std::vector<State> &path);

class Pathfinder {
    public:
        // variables
        Paths paths;
        Waiting_times mean_wt;
        double minor_pathways = 0.;
        double minor_pathways_wt = 0.;
        double forbidden_pathways = 0.;
        std::vector<State> A;             // states from
        std::vector<State> B;             // states To
        std::vector<State> C;             // forbidden states
        // set to max, to avoid conflicts in add_path for final result
        State state_from = std::numeric_limits<State>::max();
        Transition_matrix T;
        std::size_t steps = 0;                  // number of steps to propagate
        std::size_t threshold = 100000;         // number of different pathways stored
        bool isWeight = false;                  // if weight store only final state

        // debug param
        std::size_t propagated_steps = 0;
        std::size_t propagated_pathways = 0;
        double missed_final = 0.;    // did not reach final states [in percent]

        // Constructor and Destructor
        Pathfinder();
        Pathfinder(const Pathfinder& pathfinder);
        Pathfinder(const Transition_matrix& T,
                   const std::vector<State> &A,
                   const std::vector<State> &B,
                   const std::vector<State> &C,
                   const State state_from,
                   const std::size_t steps,
                   const std::size_t threshold,
                   const bool isWeight);
        virtual ~Pathfinder();

        // methods
        //! normalize the probability of all paths to 1
        void normalize();
        //! normalize the probability of all time to 1
        void normalize_times();
        //! propagate path
        virtual void propagate_path(double prob) {return;};
        //! remove unlikely paths with probability lower than cut-off
        void remove_unlikely_paths();
        //! add path to this->paths map
        void add_path(std::vector<State> path,
                      double val,
                      std::size_t threshold=10000);
        //! add path to this->paths map
        void add_time(std::vector<State> path,
                      double prob,
                      double time,
                      std::size_t threshold=10000);
        //! print pathways
        void print(int number_of_pathways=50);
        //! reset to defualt values
        void clear();
        //! return waiting time of collective pathways
        double waiting_time();
    protected:
        //! print statistics of pathways
        virtual void print_stats();
};

//! Class for running many short MCMC runs between A->B
class Mcmc : public Pathfinder {
    public:
        // variables
        std::size_t runs = 0;       // number of runs

        // Constructor and Destructor
        Mcmc();
        Mcmc(const Transition_matrix& T,
             const std::vector<State> &A,
             const std::vector<State> &B,
             const std::vector<State> &C,
             const State state_from,
             const std::size_t steps,
             const std::size_t threshold,
             const bool isWeight,
             const std::size_t runs);
        Mcmc(const Mcmc& mcmc);
        ~Mcmc();

        // methods
        //! propagate path
        void propagate_path(double prob);
    protected:
        //! print statistics of pathways
        void print_stats();
};

//! Class for running single long MCMC
class Mcmc_single : public Pathfinder {
    public:
        // variables
        std::size_t total_steps = 0;      // number of steps to propagate
        float stay_in_forbidden_region = 0.; // time of staying in forbidden region [%]

        // Constructor and Destructor
        Mcmc_single();
        Mcmc_single(const Transition_matrix& T,
                    const std::vector<State> &A,
                    const std::vector<State> &B,
                    const std::vector<State> &C,
                    const State state_from,
                    const std::size_t steps,
                    const std::size_t threshold,
                    const bool isWeight,
                    const std::size_t total_steps);
        Mcmc_single(const Mcmc_single& mcmc);
        ~Mcmc_single();

        // methods
        //! propagate path
        void propagate_path(double prob);
    protected:
        //! print statistics of pathways
        void print_stats();
};

//! Class for getting the theoretical value

class Pathfinder_tauij : public Pathfinder {
    public:
        // variables
        double cut_off = 0.;
        double error = 0.;          // percentage of not propagated pathways

        // Constructor and Destructor
        Pathfinder_tauij();
        Pathfinder_tauij(const Transition_matrix& T,
             const std::vector<State> &A,
             const std::vector<State> &B,
             const std::vector<State> &C,
             const State state_from,
             const std::size_t steps,
             const std::size_t threshold,
             const bool isWeight,
             const double cut_off = 1e-10);
        Pathfinder_tauij(const Pathfinder_tauij& pf);
        ~Pathfinder_tauij();

        // methods
        //! propagate path
        void propagate_path(double prob);
    protected:
        //! print statistics of pathways
        void print_stats();
};

#endif
