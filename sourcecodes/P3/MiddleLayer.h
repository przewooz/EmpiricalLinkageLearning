// Brian Goldman

// The Middle Layer sits between the optimization
// method and the problem that is being optimized in order to perform
// evaluation counting and progression tracking.  It records
// each time a new best fitness is reached and the evaluation at
// which it was first found.  When configured, this will also output
// the genotype of the solution which found the new best fitness.

#ifndef MIDDLELAYER_H_
#define MIDDLELAYER_H_
#include "Evaluation.h"
#include "Util.h"
#include <stdexcept>
#include <numeric>
#include <fstream>
#include "Record.h"

#include "../EATester/Log.h"

// Inherits from Evaluator in order to look like one to the optimizer
class Middle_Layer : public Evaluator {
 public:
  // Public access of recorded data
  size_t evaluations;
  float best_fitness;
  int best_found;
  std::vector<bool> best_solution;
  float best_unitation;
  Record results;
  time_t start_time;
  // Constructs a middle layer around the problem it is going to be
  // mimicing
  Middle_Layer(Configuration& config, shared_ptr<Evaluator> evaler, time_t start_time, CLog* log)
      : evaluations(0),
        best_fitness((std::numeric_limits<int>::min)()),
        best_found(0),
		best_unitation(0),
        results(config),
        outfile(config.get<string>("solution_file")),
        disable_solution_outfile(config.get<int>("disable_solution_outfile")),
        evaluator(evaler),
		start_time(start_time),
		log(log) {
  }
  // The optimization method calls this evaluation, which is a pass through
  // to the "evaluator" protected member
  virtual float evaluate(const vector<bool> & solution) override;
  float virtual max_fitness(const float max_fitness_config) override { return evaluator->max_fitness(max_fitness_config); };
  size_t virtual length(const size_t length_config) override { return evaluator->length(length_config); };

 protected:
  string outfile;
  bool disable_solution_outfile;

  shared_ptr<Evaluator> evaluator;

private:
	float unitation(const vector<bool>& solution);

	CLog* log;
};

#endif
