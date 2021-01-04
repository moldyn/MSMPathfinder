#include <iostream>
#include <iomanip>
#include <fstream>
#include "paths.hpp"
#include "main.hpp"
#include "pathfinder.hpp"
#include "transition_matrix.hpp"

#include <string>
#include <vector>
#include <cmath>

#include <chrono>

//#include <exception>
//#include <omp.h>

//! Get paths from all states a in A to B
void generate_pathfinders(Pathfinder* pathfinder,
                          Weights &weights_A){
    // generate mode specific class
    //TODO: move state_from to propagate_path function
    //TODO: remove from constructor
    // propagate pathways from all initial states
//    for(std::vector<State>::size_type i=0; i != A.size(); i++) {
    for (const auto &weight_tuple: weights_A) {
        // set initial state
        pathfinder->state_from = weight_tuple.first;

        double prob = weight_tuple.second;

        // propagete pathways
        pathfinder->propagate_path(prob);
        pathfinder->remove_unlikely_paths();
    }

    pathfinder->normalize();

    return;
}


//! Get paths from all states a in A to B
void propagate_weights(Pathfinder_map &pathfinders,
                       const std::vector<State> &A,
                       const std::vector<State> &B,
                       const std::vector<State> &C,
                       const Transition_matrix &T,
                       const double cut_off=0,
                       const std::size_t iterations=100000000,
                       const Mode mode=MCMC){
    const bool isWeight = true;
    const std::size_t steps = -1;
    const std::size_t threshold = -1;

    // set up all pathfinders to be propagated
    for(std::vector<State>::size_type i=0; i != A.size(); i++) {
        // generate mode specific class
        switch(mode){
            case PATHS:
                {
                    pathfinders[A[i]] = new Pathfinder_tauij(T, A, B, C, A[i], steps, threshold, isWeight, cut_off);
                }
                break;
            case MCMC:
                {
                    pathfinders[A[i]] = new Mcmc(T, A, B, C, A[i], steps, threshold, isWeight, iterations);
                }
                break;
            default:
                std::cerr << "    ERROR: unknown mode. this should never happen."
                          << std::endl;
                std::exit(EXIT_FAILURE);
        }
    }

    // propagate pathways and normalize them
    for (const auto pathfinder_tuple: pathfinders) {
        pathfinder_tuple.second->propagate_path(1.);
        pathfinder_tuple.second->normalize();
    }
    return;
}


//! MAIN ROUTINE
int run_paths(Mode mode,
              std::vector<State> A,
              std::vector<State> B,
              std::vector<State> C,
              double cut_off,
              std::size_t steps,
              std::size_t iterations,
              std::size_t total_steps,
              Transition_matrix T,
              std::string output,
              std::size_t threshold,
              int argc, char* argv[]){


    Pathfinder_map pathfinders_from_A_weights, pathfinders_from_B_weights;
    // needed to calculate weights
    // MCMC mode because of better initial convergence
    if ((A.size() > 1) && (mode != MCMCSINGLE)) {
        // remove diagonal
        bool keep_diag = T.get_keep_diag();
        T.set_keep_diag(false);

        std::cout << "\n~~~ Get paths for N->infty for weighting" << std::endl;
        if (B.size() > 1) {
            propagate_weights(pathfinders_from_A_weights, A, B, C, T);
        }
        propagate_weights(pathfinders_from_B_weights, B, A, C, T);

        // restore diagonal setting
        T.set_keep_diag(keep_diag);
    }

    // calculate weights
    // TODO: use matrix form
    Weights weights_A;
    if ((A.size() > 1) && (B.size() == 1) && (mode != MCMCSINGLE)) {
        // initialize weights with a zeros
        for (const auto &state: A) {
            weights_A[state] = 0.;
        }
        // get weights
        for (const auto &pathfinder_tuple: pathfinders_from_B_weights) {
            for (const auto &path_tuple: pathfinder_tuple.second->paths) {
                weights_A[path_tuple.first.back()] += path_tuple.second;
            }
        }
    } else if ((A.size() > 1) && (B.size() > 1) && (mode != MCMCSINGLE)) {
        std::cout << "\n" << std::endl;
        // initialize weights equally
        Weights weights_B, prev_weights_A;
        for (const auto &state: A) {
            weights_A[state] = 1./A.size();
        }
        for (const auto &state: B) {
            weights_B[state] = 1./B.size();
        }

        double diff = 1;
        double diff_old = diff;
        // propagate from A->B->A and check how much weight changed.
        int iteration_converged = 0;
        Progressspinner spinner("~~~ Calculate weights");
        int i=0;
        for(; i<300; ++i) {
            spinner.next();
            for (const auto &state: A) {
                prev_weights_A[state] = weights_A[state];
            }

            // set weights_B to zero
            for (const auto &state_weight: weights_B){
                weights_B[state_weight.first] = 0.;
            }

            // find weights_B for probability weights_A to start in A
            for (const auto &state_weight: weights_A){
                // state stateWeight.first
                // weight stateWeight.second
                for (const auto &pathTuple: pathfinders_from_A_weights[state_weight.first]->paths){
                    weights_B[pathTuple.first.back()] +=
                        pathTuple.second*state_weight.second;
                }
            }

            // set weights_A to zero
            for (const auto &state_weight: weights_A){
                weights_A[state_weight.first] = 0.;
            }


            // find weights_A for probability weights_B to start in B
            for (const auto &state_weight: weights_B){
                // state stateWeight.first
                // weight stateWeight.second
                for (const auto &pathTuple: pathfinders_from_B_weights[state_weight.first]->paths){
                    weights_A[pathTuple.first.back()] +=
                        pathTuple.second*state_weight.second;
                }
            }

            // normalize weight
            double sum = 0.;
            for (const auto &state: A) {
                sum += weights_A[state];
            }
            for (const auto &state: A) {
                weights_A[state] /= sum;
            }

            diff_old = diff;
            diff = 0.;
            for (const auto &state: A) {
                diff += std::fabs(weights_A[state] - prev_weights_A[state]);
            }

            if (diff_old == diff) {
                ++iteration_converged;
                if (iteration_converged > 25) {
                    break;
                }
            }

            if (diff < 1e-20) {
                break;
            }
        }
        spinner.close();
        std::cout << "    iterations to find weights: " << i
                  << " with error " << diff << std::endl;

    } else {
        weights_A[A[0]] = 1.;
    }


    // loop over all a in A
    std::cout << "~~~ Get paths from initial states" << std::endl;

    // start timer
    std::chrono::steady_clock::time_point t_start, t_end;
    t_start = std::chrono::steady_clock::now();

    Pathfinder* pathfinder_from_A;
    const bool isWeight = false;
    switch(mode){
        case PATHS:
            {
                pathfinder_from_A = new Pathfinder_tauij(T, A, B, C, A[0], steps, threshold, isWeight, cut_off);
            }
            break;
        case MCMC:
            {
                pathfinder_from_A = new Mcmc(T, A, B, C, A[0], steps, threshold, isWeight, iterations);
            }
            break;
        case MCMCSINGLE:
            {
                pathfinder_from_A = new Mcmc_single(T, A, B, C, A[0], steps, threshold, isWeight, total_steps);
            }
            break;
        default:
            std::cerr << "    ERROR: unknown mode. this should never happen."
                      << std::endl;
            std::exit(EXIT_FAILURE);
    }

    // generate mode specific class
    generate_pathfinders(pathfinder_from_A,
                         weights_A);

    std::cout<< "\n\n~~~ FINAL PATHWAYS:" << std::endl;
    pathfinder_from_A->print();

    t_end = std::chrono::steady_clock::now();
    float t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count()/1000.;


    //TODO: use same function for print and storing
    if( output[0] != '\0' ) {  // check if string is non empty
        // sort
        Paths_sorted final_paths_sorted;
        for (const auto &path_tuple: pathfinder_from_A->paths) {
            final_paths_sorted.insert({path_tuple.second, path_tuple.first});
        }

        std::cout << "\n\n~~~ save output" << std::endl;
        std::ofstream output_file (output);
        output_file << std::fixed
                    << std::setprecision(10);
        if (output_file.is_open()) {
            Paths_sorted::reverse_iterator it_r;
        	it_r = final_paths_sorted.rbegin();
        	std::size_t i = 0;
        	std::size_t i_max = threshold;
            double count = 0;
            int first = -1;
            // header comment
        	output_file << "# This file was created by MSMPathfinder with:\n#\n# ";
            std::vector<std::string> arguments(argv, argv + argc);
            for (std::string& arg_string : arguments){
                output_file << arg_string << " ";
            }
            output_file << "\n";
            switch(mode){
                case MCMC:
                    break;
                case MCMCSINGLE:
                    break;
                case PATHS:
                    output_file << "# undiscovered pathspace [%]: "
                                << 100*static_cast<Pathfinder_tauij*>(pathfinder_from_A)->error
                                << "\n";
                    break;
            }
  	        output_file << "# minor pathways [%]: "
                        << 100*pathfinder_from_A->minor_pathways << "\n"
                        << "# minor pathways wt [t_lag]: "
                        << pathfinder_from_A->minor_pathways_wt << "\n"
                        << "# missed pathways [%]: "
                        << pathfinder_from_A->missed_final << "\n"
                        << "# number of propagated steps: "
                        << pathfinder_from_A->propagated_steps << "\n"
                        << "# propagated pathways: "
                        << pathfinder_from_A->propagated_pathways << "\n"
                        << "# waiting time [t_lag]: "
                        << pathfinder_from_A->waiting_time() << "\n"
                        << "# number of different pathways: "
                        << pathfinder_from_A->paths.size() << "\n"
                        << "# elapsed time [s]: "
                        << t_elapsed << "\n"
        	            << "# count [%]  total [%] time [t_lag] pathway\n";
            for (; it_r != final_paths_sorted.rend(); it_r++){

                if (first == -1) {
                    first = it_r->first;
                } else if (i >= i_max) {
                	output_file << " # ...  ...  ...  ...  ...\n";
                    break;
                }
                count += 100*it_r->first;
                output_file << 100*it_r->first << " "
                            << count << " "
                            << pathfinder_from_A->mean_wt[it_r->second] << " ";
                for (const auto& state : it_r->second) {
                    output_file << state+1 << " ";
                }
                output_file << "\n";
                ++i;
            }
            output_file.close();
        } else {
            std::cout << "    ERROR: Unable to open file" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    //TODO: delete all objects created with new
    std::exit(EXIT_SUCCESS);
}
