#include <iostream>
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
void generate_pathfinders(const Mode mode,
                          Pathfinder_map &pathfinders,
                          const std::vector<std::size_t> &A,
                          const std::vector<std::size_t> &B,
                          const std::vector<std::size_t> &C,
                          const Transition_matrix &T,
                          const std::size_t steps,
                          const std::size_t total_steps,
                          const std::size_t iterations,
                          const double cut_off,
                          const std::size_t threshold,
                          const bool isWeight,
                          const bool normalize){
    // set up all pathfinders to be propagated
    for(std::vector<std::size_t>::size_type i=0; i != A.size(); i++) {
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
            case MCMCSINGLE:
                if (i==0) {
                    pathfinders[A[i]] = new Mcmc_single(T, A, B, C, A[i], steps, threshold, isWeight, total_steps);
                }
                break;
            default:
                std::cerr << "    ERROR: unknown mode. this should never happen."
                          << std::endl;
                std::exit(EXIT_FAILURE);
        }
    }

    // propagate pathways
    for (const auto pathfinder_tuple: pathfinders) {
        pathfinder_tuple.second->propagate_path();
//        std::cout << "                 unique pathways: "
//                  << pathfinder_tuple.second->paths.size();
        pathfinder_tuple.second->remove_unlikely_paths();
//        std::cout << "/" << pathfinder_tuple.second->paths.size()
//                  << " with "
//                  << 100.*(1.-pathfinder_tuple.second->minor_pathways)
//                  << "%" <<std::endl;
        if (normalize) {
            pathfinder_tuple.second->normalize();
        }
    }
    return;
}


//! MAIN ROUTINE
int run_paths(Mode mode,
              std::vector<std::size_t> A,
              std::vector<std::size_t> B,
              std::vector<std::size_t> C,
              double cut_off,
              std::size_t steps,
              std::size_t iterations,
              std::size_t total_steps,
              Transition_matrix T,
              std::string output,
              std::size_t threshold,
              int argc, char* argv[]){

    Pathfinder_map pathfinders_from_A;
    Pathfinder_map pathfinders_from_A_weights, pathfinders_from_B_weights;

    // loop over all a in A
    std::cout << "~~~ Get paths from initial states" << std::endl;

    // start timer
    std::chrono::steady_clock::time_point t_start, t_end;
    t_start = std::chrono::steady_clock::now();

    generate_pathfinders(mode,
                         pathfinders_from_A,
                         A, B, C, T,
                         steps,
                         total_steps,
                         iterations,
                         cut_off,
                         threshold,
                         false,
                         false);

    // TODO:DEBUG
//    for(const auto &pathfinder_tuple: pathfinders_from_A) {
//        std::cout << "~~~ Paths from state: " << pathfinder_tuple.first << std::endl;
//        pathfinder_tuple.second->print();
//    }
    // END DEBUG

    // needed to calculate weights
    if ((A.size() > 1) && (mode != MCMCSINGLE)) {
        // remove diagonal
        T.set_keep_diag(false);
        std::cout << "\n~~~ Get paths for N->infty for weighting" << std::endl;
        if (B.size() > 1) {
            // MCMC mode because of better initial convergence
           generate_pathfinders(MCMC,
                                pathfinders_from_A_weights,
                                A, B, C, T,
                                -1,  // fixed number of steps
                                0,//total_steps/10,
                                100000000,  //iterations/10,
                                0,  //cut_off*10,
                                threshold,
                                true,
                                true);
             // generate_pathfinders(PATHS,
             //                      pathfinders_from_A_weights,
             //                      A, B, C, T,
             //                      -1,  // fixed number of steps
             //                      0,//total_steps/10,
             //                      0,  //iterations/10,
             //                      0.00000000001,  //cut_off*10,
             //                      threshold,
             //                      true,
             //                      true);
        }
        // MCMC mode because of better initial convergence
        generate_pathfinders(MCMC,
                             pathfinders_from_B_weights,
                             B, A, C, T,
                             -1,  // fixed number of steps
                             0,//total_steps/10,
                             100000000,  //iterations/10,
                             0,  //cut_off*10,
                             threshold,
                             true,
                             true);
          // generate_pathfinders(PATHS,
          //                      pathfinders_from_B_weights,
          //                      B, A, C, T,
          //                      -1,  // fixed number of steps
          //                      0,//total_steps/10,
          //                      0,  //iterations/10,
          //                      0.00000000001,  //cut_off*10,
          //                      threshold,
          //                      true,
          //                      true);
        // TODO:DEBUG
//        for(const auto &pathfinder_tuple: pathfinders_from_B) {
//            std::cout << "~~~ Paths from state: " << pathfinder_tuple.first << std::endl;
//            pathfinder_tuple.second->print();
//        }
        // END DEBUG
    }

    // calculate weights
    Weights weights_A;
    if ((A.size() > 1) && (B.size() == 1) && (mode != MCMCSINGLE)) {
        // initialize weights with a zeros
        for (const auto& state: A) {
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
        for (const auto& state: A) {
            weights_A[state] = 1./A.size();
        }
        for (const auto& state: B) {
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
                weights_B[state_weight.first] = 0;
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
                weights_A[state_weight.first] = 0;
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
            double sum = 0;
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
        for (const auto &pathfinder_tuple: pathfinders_from_A) {
            weights_A[pathfinder_tuple.first] = 1;
        }
    }

    // free memory of weight calculations
//    for (auto pathfinder_tuple: pathfinders_from_A_weights) {
//        delete pathfinder_tuple.second;
//    }
//    for (auto pathfinder_tuple: pathfinders_from_B_weights) {
//        delete pathfinder_tuple.second;
//    }

    // merge weights into final paths
    Pathfinder* paths_result;
    switch(mode){
        case PATHS:
            {
                paths_result = new Pathfinder_tauij();
            }
            break;
        case MCMC:
            {
                paths_result = new Mcmc();
            }
            break;
        case MCMCSINGLE:
            {
                paths_result = new Mcmc_single();
            }
            break;
        default:
            std::cerr << "    ERROR: unknown mode. this should never happen."
                      << std::endl;
            return EXIT_FAILURE;
    }

    for (const auto &pathfinder_tuple: pathfinders_from_A) {
        for (const auto &path_tuple: pathfinder_tuple.second->paths) {
            paths_result->add_path(path_tuple.first,
                                   path_tuple.second*weights_A[pathfinder_tuple.first],
                                   threshold);
            paths_result->add_time(path_tuple.first,
                                   path_tuple.second*weights_A[pathfinder_tuple.first],
                                   pathfinder_tuple.second->mean_wt[path_tuple.first],
                                   threshold);
        }
        paths_result->propagated_steps +=
            pathfinder_tuple.second->propagated_steps;
        paths_result->minor_pathways +=
            pathfinder_tuple.second->minor_pathways*weights_A[pathfinder_tuple.first];
        paths_result->missed_final +=
            pathfinder_tuple.second->missed_final*weights_A[pathfinder_tuple.first];
        paths_result->propagated_pathways +=
            pathfinder_tuple.second->propagated_pathways;
        switch(mode){
            case PATHS:
                static_cast<Pathfinder_tauij*>(paths_result)->error +=
                    static_cast<Pathfinder_tauij*>(pathfinder_tuple.second)->error
                        *weights_A[pathfinder_tuple.first];
                break;
            case MCMC:
                break;
            case MCMCSINGLE:
                static_cast<Mcmc_single*>(paths_result)->stay_in_forbidden_region +=
                    static_cast<Mcmc_single*>(pathfinder_tuple.second)->stay_in_forbidden_region;
                break;
            default:
                std::cerr << "    ERROR: unknown mode. this should never happen."
                          << std::endl;
                return EXIT_FAILURE;
        }
    }

    // free memory of weight calculations
//    for (auto pathfinder_tuple: pathfinders_from_A) {
//        delete pathfinder_tuple.second;
//    }

    std::cout<< "\n\n~~~ FINAL PATHWAYS:" << std::endl;
    paths_result->normalize_times();
    paths_result->normalize();
    paths_result->print();

    //TODO: end timer
    t_end = std::chrono::steady_clock::now();
    float t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count()/1000.;


    // TODO: store results
    if( output[0] != '\0' ) {  // check if string is non empty
        // sort
        Paths_sorted final_paths_sorted;
        for (const auto &path_tuple: paths_result->paths) {
            final_paths_sorted.insert({path_tuple.second, path_tuple.first});
        }

        std::cout << "\n\n~~~ save output" << std::endl;
        std::ofstream output_file (output);
        output_file << std::fixed;
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
                                << 100*static_cast<Pathfinder_tauij*>(paths_result)->error
                                << "\n";
                    break;
            }
  	        output_file << "# minor pathways [%]: "
                        << 100*paths_result->minor_pathways << "\n"
                        << "# missed pathways [%]: "
                        << paths_result->missed_final << "\n"
                        << "# number of propagated steps: "
                        << paths_result->propagated_steps << "\n"
                        << "# propagated pathways: "
                        << paths_result->propagated_pathways << "\n"
                        << "# number of different pathways: "
                        << paths_result->paths.size() << "\n"
                        << "# elapsed time [s]: "
                        << t_elapsed << "\n"
        	            << "# count [%]  total [%] time [steps] pathway\n";
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
                            << paths_result->mean_wt[it_r->second] << " ";
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