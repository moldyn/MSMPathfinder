#include <iostream>
#include "main.hpp"
#include "pathfinder.hpp"

#include <cmath>
#include <algorithm>
#include <random>
// Default constructor
Pathfinder::Pathfinder() {}
Mcmc::Mcmc() {}
Mcmc_single::Mcmc_single() {}
Pathfinder_tauij::Pathfinder_tauij() {}

// Overloaded constructor
Pathfinder::Pathfinder(const Transition_matrix& T,
                       const std::vector<State> &A,
                       const std::vector<State> &B,
                       const std::vector<State> &C,
                       const State state_from,
                       const std::size_t steps,
                       const std::size_t threshold,
                       const bool isWeight)
    : minor_pathways(0.), minor_pathways_wt(0.), forbidden_pathways(0.),
      A(A), B(B), C(C), state_from(state_from), T(T), steps(steps),
      threshold(threshold), isWeight(isWeight), propagated_steps(0.),
      propagated_pathways(0), missed_final(0.)
{}
Mcmc::Mcmc(const Transition_matrix& T,
           const std::vector<State> &A,
           const std::vector<State> &B,
           const std::vector<State> &C,
           const State state_from,
           const std::size_t steps,
           const std::size_t threshold,
           const bool isWeight,
           const std::size_t runs)
    : Pathfinder(T, A, B, C, state_from, steps, threshold, isWeight),
      runs(runs)
{}
Mcmc_single::Mcmc_single(const Transition_matrix& T,
                         const std::vector<State> &A,
                         const std::vector<State> &B,
                         const std::vector<State> &C,
                         const State state_from,
                         const std::size_t steps,
                         const std::size_t threshold,
                         const bool isWeight,
                         const std::size_t total_steps)
    : Pathfinder(T, A, B, C, state_from, steps, threshold, isWeight),
      total_steps(total_steps), stay_in_forbidden_region(0.)
{}
Pathfinder_tauij::Pathfinder_tauij(const Transition_matrix& T,
                                   const std::vector<State> &A,
                                   const std::vector<State> &B,
                                   const std::vector<State> &C,
                                   const State state_from,
                                   const std::size_t steps,
                                   const std::size_t threshold,
                                   const bool isWeight,
                                   const double cut_off)
    : Pathfinder(T, A, B, C, state_from, steps, threshold, isWeight),
      cut_off(cut_off), error(0.)
{}

// copy constructor
Pathfinder::Pathfinder(const Pathfinder& pathfinder)
    : paths(pathfinder.paths),
      mean_wt(pathfinder.mean_wt),
      minor_pathways(pathfinder.minor_pathways),
      minor_pathways_wt(pathfinder.minor_pathways_wt),
      forbidden_pathways(pathfinder.forbidden_pathways),
      A(pathfinder.A),
      B(pathfinder.B),
      C(pathfinder.C),
      state_from(pathfinder.state_from),
      T(pathfinder.T),
      steps(pathfinder.steps),
      threshold(pathfinder.threshold),
      isWeight(pathfinder.isWeight),
      propagated_steps(pathfinder.propagated_steps),
      propagated_pathways(pathfinder.propagated_pathways),
      missed_final(pathfinder.missed_final)
{}
Mcmc::Mcmc(const Mcmc& mcmc)
    : Pathfinder(mcmc), runs(mcmc.runs)
{}
Mcmc_single::Mcmc_single(const Mcmc_single& mcmc)
    : Pathfinder(mcmc), total_steps(mcmc.total_steps),
      stay_in_forbidden_region(mcmc.stay_in_forbidden_region)
{}
Pathfinder_tauij::Pathfinder_tauij(const Pathfinder_tauij& pf)
    : Pathfinder(pf), cut_off(pf.cut_off), error(pf.error)
{}

// Destructor
Pathfinder::~Pathfinder() {this->T.~Transition_matrix();}
Mcmc::~Mcmc() {}
Mcmc_single::~Mcmc_single() {}
Pathfinder_tauij::~Pathfinder_tauij() {}

// Propagate Paths
void Mcmc::propagate_path(double prob=1.)
{
    // allocate memory for path initialize it
    std::vector<State> path;
    if (this->steps > 1000) {  // maximal length of reserved storage
        path.reserve(1000);
    } else {
        path.reserve(this->steps);
    }
    double rand = 0.;

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1.0);

    Progressbar eta_bar(this->runs, this->state_from+1);
    for (std::size_t run=0; run<this->runs; ++run) {
        path.clear();  // reset path
        State next_state = this->state_from;
        path.push_back(next_state);
        for(std::size_t step=0; step<this->steps; ++step) {
            rand = dis(gen);
            next_state = this->T.get_mcmc(next_state, rand);
            path.push_back(next_state);
            // check if loop and remove
            simplify_path(path);

            this->propagated_steps++;

            // check if reached forbidden state
            if (std::count(this->C.begin(), this->C.end(), next_state)) {
                this->forbidden_pathways += prob/this->runs;
                goto reached_final_state;
            }
            // check if reached final state
            if (std::count(this->B.begin(), this->B.end(), next_state)) {
                this->add_path(path, prob/this->runs, -1);  //this->threshold);
                this->add_time(path, prob/this->runs, step + 1, -1); //this->threshold);
                this->propagated_pathways++;
                goto reached_final_state;
            }
            // check if returned to starting state
            if (std::count(this->A.begin(), this->A.end(), next_state)) {
                path.clear();
                path.push_back(next_state);
            }
        }

        // did not reach one of the final states B
        this->missed_final += prob/this->runs; ///this->runs;

        // print missed pathways
//        for (const auto& state : path) {
//            std::cout << state+1 << " ";
//        }
//        std::cout << "\b \b"  // remove last ->
//                  << std::endl;

        reached_final_state: ; // goto for breaking nested loop
        eta_bar.next();
    }
    eta_bar.close();

    // normalize missed final to percentage
    this->missed_final *= prob*100./this->runs;

    // normalize times
    //this->normalize_times();

    return;
}

// Propagate Paths
void Mcmc_single::propagate_path(double prob=1.)
{
    // allocate memory for path initialize it
    std::vector<State> path;
    if (this->steps > 1000) {  // maximal length of reserved storage
        path.reserve(1000);
    } else {
        path.reserve(this->steps);
    }
    double rand = 0.;
    path.clear();  // reset path

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1.0);

    bool propagating_forward = true;  // always starting in A

    State next_state = this->state_from;
    path.push_back(next_state);
    std::size_t step_started = 0;

    Progressbar eta_bar(this->total_steps);
    for(std::size_t step=0; step<this->total_steps; ++step) {
        rand = dis(gen);
        next_state = this->T.get_mcmc(next_state, rand);
        if ((path.size() == 0) || (next_state != path.back())) {
            path.push_back(next_state);
        }
        // check if loop and remove
        simplify_path(path);

        // check if reached forbidden state
        if (std::count(this->C.begin(), this->C.end(), next_state)) {
            this->stay_in_forbidden_region += 1.;
            path.clear();
            propagating_forward = false;
        }
        // check if reached final state
        if (std::count(this->B.begin(), this->B.end(), next_state)) {
            // if path is too long, throw it as well
            if (propagating_forward) {
                if (step-step_started <= this->steps) {
                    this->add_path(path, 1., this->threshold);
                    this->add_time(path, 1., step-step_started, this->threshold); // ;-1);
                    this->propagated_pathways++;
                    propagating_forward = false;
                } else {
                    // did not reach within number of given steps
                    path.clear();
                    propagating_forward = false;
                    this->missed_final += 1.;
                }
            }
        }
        // check if returned to starting state
        if (std::count(this->A.begin(), this->A.end(), next_state)) {
            path.clear();
            path.push_back(next_state);
            // this is for counting from first entering the region A
            if (!propagating_forward) {
                step_started = step;
                propagating_forward = true;
            }
        }
        eta_bar.next();
    }
    eta_bar.close();

    this->propagated_steps = this->total_steps;
    // normalize missed final to percentage
    this->stay_in_forbidden_region *= 100./this->total_steps;
    this->missed_final *= 100./(this->propagated_pathways + this->missed_final);

    // normalize times
    // this->normalize_times();

    return;
}

void Pathfinder_tauij::propagate_path(double prob=1.)
{
   // allocate memory for path initialize it
    std::vector<State> path;
    if (this->steps > 1000) {  // maximal length of reserved storage
        path.reserve(1000);
    } else {
        path.reserve(this->steps);
    }
    Pathfinder pathfinder_curr, pathfinder_prev;

    State current_state=this->state_from;
    path.push_back(this->state_from);

    Progressspinner spinner(this->state_from+1, 100000);
    //Progressbar eta_bar(this->runs, this->state_from);
    pathfinder_prev.add_path(path, prob, this->threshold);

    for(std::size_t step=0; step<this->steps; ++step) {
        pathfinder_curr.clear();  // empty all
        if (pathfinder_prev.paths.empty()){
            break;
        }
        for(const auto paths_tuple : pathfinder_prev.paths) {
            spinner.next();

            current_state = paths_tuple.first.back();
            // append all possible paths
            for(std::size_t state=0; state<this->T.get_no_of_states(); state++) {
                prob = paths_tuple.second*this->T.get_Tij(current_state, state);

                if (prob > this->cut_off) {
                    // check if reached forbidden state
                    if (std::count(this->C.begin(), this->C.end(), state)) {
                        this->forbidden_pathways += prob;
                    } else {
                        // check if returned to starting state
                        if (std::count(this->A.begin(), this->A.end(), state)) {
                         path.clear();
                        } else {
                         path = paths_tuple.first;
                        }
                        // DEBUG
                        //path = paths_tuple.first;
                        // DEBUG

                        if ((path.size() == 0) || (state != path.back())) {
                            path.push_back(state);
                        }
                        simplify_path(path);
                        // check if reached final state
                        if (std::count(this->B.begin(), this->B.end(), state)) {
                            this->add_path(path, prob, -1); //this->threshold);
                            this->add_time(path, prob, step+1, -1); //this->threshold);
                            this->propagated_pathways++;
//                        } else if (pathfinder_curr.paths.size() > this->threshold) {
//                            // add to minor pathways
//                            this->minor_pathways += prob;
//                            this->minor_pathways_wt += prob*(step+1);
                        } else {
                            pathfinder_curr.add_path(path, prob, -1); //this->threshold);
                        }
                    }
                } else {
                    this->error += prob;
                }
            }
        }
        pathfinder_prev = pathfinder_curr;
    }

    for(const auto paths_tuple : pathfinder_prev.paths) {
        this->missed_final += 100*paths_tuple.second;
    }

    // free memory
    pathfinder_curr.clear();
    pathfinder_prev.clear();

    // normalize times
    //this->normalize_times();

    spinner.close();

    return;
}


double Pathfinder::waiting_time()
{
    double wt = 0.;
    for (auto& path : this->mean_wt) {
        wt += path.second*this->paths[path.first];
    }
    wt += this->minor_pathways_wt*this->minor_pathways;
    return wt;
}

// TODO: add normalization to prob!=1.
void Pathfinder::normalize_times()
{
    // normalize waiting times
    for (auto& path : this->mean_wt) {
        path.second /= this->paths[path.first];
    }
    if (this->minor_pathways > 0) {
        this->minor_pathways_wt /= this->minor_pathways;
    }
    return;
}


// normalize paths
// TODO: add normalization to prob!=1.
void Pathfinder::normalize()
{
    double sum = 0;
    for (const auto& path : this->paths) {
        sum += path.second;
    }
    sum += this->minor_pathways;

    if (sum > 0){  // if at least a single entry normalize
        for (auto& path : this->paths) {
            path.second /= sum;
        }
        this->minor_pathways /= sum;
        this->forbidden_pathways /= sum;

        // normalize times
        for (auto& path : this->mean_wt) {
            path.second /= sum * this->paths[path.first];
        }
        if (this->minor_pathways > 0) {
            this->minor_pathways_wt /= sum * this->minor_pathways;
        }
    }
    return;
}


// removes loop between last state and some intermediate
void simplify_path(std::vector<State>& path)
{
    if ((path.empty()) || (path.size() == 1)) {
        return;
    }
    // Check if last element exists already in path
    auto it = std::find(path.begin(), --path.end(), path.back());

    if (it == --path.end()) {
        return;
    } else {
        // remove loop
        int index = std::distance(path.begin(), it);
        path.erase(path.begin()+index, --path.end());
    }
    return;
}

//! add path to paths map
void Pathfinder::add_path(std::vector<State> path,
                          double val,
                          std::size_t threshold)
{
// if returned add to self transition probability
//    if (path.back() == this->state_from) {
//        std::cerr << " eventually add backtransition to self transition rate"
//                  << std::endl;
//    }
    if (isWeight) {
        path = {path.back()};
    }
    Paths::iterator it(this->paths.find(path));
    if (it != this->paths.end()){
        it->second += val;
    } else {
        if (this->paths.size() <= threshold) {
            this->paths[path] = val;
        } else {
            this->minor_pathways += val;
        }
    }
    return;
}

//! add path to paths map
void Pathfinder::add_time(std::vector<State> path,
                          double prob,
                          double time,
                          std::size_t threshold)
{
    if (isWeight) {
        path = {path.back()};
    }
    Paths::iterator it(this->mean_wt.find(path));
    if (it != this->mean_wt.end()){
        it->second += prob*time;
    } else {
        if (this->mean_wt.size() <= threshold) {
            this->mean_wt[path] = prob*time;
        } else {
            this->minor_pathways_wt += prob*time;
        }
    }
    return;
}


// remove most unlike paths
void Pathfinder::remove_unlikely_paths()
{

    // keep only first nth paths
    if (this->paths.size() > this->threshold) {
        // sort pathways
        Paths_sorted final_paths_sorted;
        for (const auto &path_tuple: this->paths) {
            final_paths_sorted.insert({path_tuple.second, path_tuple.first});
        }

        // generate iterater decreasingly ordered
        Paths_sorted::reverse_iterator it_r;
        it_r = final_paths_sorted.rbegin();

        // clear paths
        this->paths.clear();

        for(std::size_t i=0; it_r != final_paths_sorted.rend(); ++i) {
            if(i >= this->threshold) {
                this->minor_pathways += it_r->first;
            } else {
                this->paths[it_r->second] = it_r->first;
            }
            ++it_r;
        }

//        Paths::iterator it = this->paths.begin();
//        Paths::iterator endIt = this->paths.end();
//        for(std::size_t i=0; it != endIt; ++i) {
//            if(i > this->threshold) {
//                this->minor_pathways += it->second;
//                this->paths.erase(it++);
//            } else {
//                ++it;
//            }
//        }
    }
    // normalize result
    // this->normalize();
    return;
}


void Pathfinder::print(int number_of_pathways)
{
    Paths_sorted final_paths_sorted;
    for (const auto &path_tuple: this->paths) {
        final_paths_sorted.insert({path_tuple.second, path_tuple.first});
    }

    // print class sepcific statistics
    this->print_stats();
    // print general statistics
    std::cout << "    missed pathways [%]: "
              << this->missed_final << "\n"
              << "    minor pathways [%]: "
              << this->minor_pathways*100 << "\n"
              << "    minor pathways wt [t_lag]: "
              << this->minor_pathways_wt << "\n"
              << "    waiting time [t_lag]: "
              << this->waiting_time() << "\n"
              << "    propagated pathways: "
              << this->propagated_pathways << "\n"
              << "    number of different pathways: "
              << this->paths.size() << "\n" << std::endl;

    // print result
    Paths_sorted::reverse_iterator it_r;
    it_r = final_paths_sorted.rbegin();
    int j = 0;
    double count = 0;
    double first = -1;
    std::cout << "    count      of 1st     total      wt [t_lag]  pathway" << std::endl;
    std::cout << std::fixed;
    for (; it_r != final_paths_sorted.rend(); it_r++){
        if (first == -1) {
            first = it_r->first;
        }
        count += it_r->first;
        std::cout << "    " << it_r->first
                  << "   " << it_r->first/first
                  << "   " << count
                  << "   " << this->mean_wt[it_r->second] << "   ";
        for (const auto& state : it_r->second) {
            std::cout << state+1 << "->";
        }
        std::cout << "\b\b  \b\b"  // remove last ->
                  << std::endl;
        ++j;
        if (j >= number_of_pathways){
        	std::cout << "    ...        ...        ...        ..." << std::endl;
            break;
        }
    }
    std::cout << std::endl;
    return;
}

void Pathfinder::print_stats()
{
    // total counts are normalized to error+1, so total error is...
    std::cout << "    forbidden pathways: "
              << this->forbidden_pathways  << "\n";
    return;
}

void Mcmc::print_stats()
{
    // total counts are normalized to error+1, so total error is...
    std::cout << "    forbidden pathways: "
              << this->forbidden_pathways  << "\n";
    return;
}

void Mcmc_single::print_stats()
{
    // total counts are normalized to error+1, so total error is...
    std::cout << "    number of propagated steps: "
              << this->propagated_steps << "\n"
              << "    stay in forbidden region [%]: "
              << this->stay_in_forbidden_region  << "\n";
    return;
}

void Pathfinder_tauij::print_stats()
{
    // total counts are normalized to error+1, so total error is...
    std::cout << "    with missing pathspace [%]: "
              << 100*this->error << "\n"
              << "    forbidden pathways: "
              << this->forbidden_pathways  << "\n";
    return;
}

void Pathfinder::clear()
{
    this->paths.clear();
    this->A.clear();
    this->B.clear();
    this->C.clear();
    this->T.clear();
    this->isWeight = false;
    this->minor_pathways = 0.;
    this->forbidden_pathways = 0.;
    this->state_from = std::numeric_limits<State>::max();
    this->propagated_steps = 0;
    this->propagated_pathways = 0;
    this->missed_final = 0.;
    return;
}
