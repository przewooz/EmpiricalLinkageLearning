/*
 * Experiments.cpp
 *
 *  Created on: Jul 20, 2013
 *      Author: goldman
 */

#include "Experiments.h"

// A single run of optimization
Record single_run(Random& rand, Configuration& config,
                  evaluation::pointer problem, optimize::pointer solver,
                  int run, CLog* log) {
  time_t time_limit = config.get<int>("time_limit");
  time_t start_time = time(nullptr);

  // Middle Layer's sit between the problem and the solver, tracking optimization
  Middle_Layer recorder(config, problem(config, run), start_time, log);
  auto optimizer = solver(rand, recorder, config);

  float good_enough = recorder.max_fitness(config.get<int>("fitness_limit"));

  // Iterate the optimizer until the solution is reached, the maximum number
  // of evaluations is performed, or the optimizer reaches stagnation
  bool improvement_possible = true;

  time_t elapsed_time;

  do {
	improvement_possible = optimizer->iterate();
	elapsed_time = time(nullptr) - start_time;
  } while (recorder.best_fitness < good_enough && elapsed_time < time_limit
	  && improvement_possible);

  return recorder.results;
}

// Performs multiple runs of optization
vector<Record> multirun(Random& rand, Configuration& config,
                        evaluation::pointer problem, optimize::pointer solver) {
  int runs = config.get<int>("runs");
  int verbosity = config.get<int>("verbosity");
  vector<Record> records;
  // perform the desired number of runs
  for (int run = 0; run < runs; run++) {
    records.push_back(single_run(rand, config, problem, solver, run, nullptr));
    // output partial summaries
    if (verbosity > 0) {
      auto summary = Record::summarize(records, config);
      std::cout << "Run: " << run << " Evals: "
                << std::get<1>(records[records.size() - 1].best()) << " MES: "
                << summary[MES] << " MAD: " << summary[MAD] << " FAILURES: "
                << summary[FAILURES] << std::endl;
    }
  }
  return records;
}

// Determine the minimum useful population size
int bisection(Random& rand, Configuration& config, evaluation::pointer problem,
              optimize::pointer solver) {
  int runs = config.get<int>("runs");
  float good_enough = config.get<int>("fitness_limit");
  vector<Record> records;
  // initial bounds
  int min = 0;
  int max = 1;
  size_t failed_on = 0;
  bool success;
  // Double the maximum size until a successful run is found
  do {
    // double the bounds each time
    min = max;
    max *= 2;
    config.set("pop_size", max);
    std::cout << "Pop size: " << max << std::endl;
    success = true;
    // Perform multiple runs, stopping as soon as one of the runs failed to find
    // the global optimum
    for (int i = 0; success && i < runs; i++) {
      int problem_number = ((i + failed_on) % runs);
      std::cout << "\tTrying problem: " << problem_number << std::endl;
      auto results = single_run(rand, config, problem, solver, problem_number, nullptr);
      if (std::get<0>(results.best()) < good_enough) {
        failed_on = problem_number;
        success = false;
      }
    }
  } while (!success);

  // Bisect between the min (unsuccessful) and max (successful) until they meet
  int guess;
  while (min + 1 < max) {
    guess = (max + min) / 2;
    std::cout << "Pop size: " << guess << std::endl;
    config.set("pop_size", guess);
    success = true;
    // Perform multiple runs, stopping as soon as one of the runs failed to find
    // the global optimum
    for (int i = 0; success && i < runs; i++) {
      int problem_number = ((i + failed_on) % runs);
      std::cout << "\tTrying problem: " << problem_number << std::endl;
      auto results = single_run(rand, config, problem, solver, problem_number, nullptr);
      if (std::get<0>(results.best()) < good_enough) {
        failed_on = problem_number;
        success = false;
      }
    }
    // update the bound
    if (success) {
      max = guess;
    } else {
      min = guess;
    }
  }
  return max;
}
