// Brian Goldman

// Implements a wrapper around an evaluation function to
// allow for fitness trajectory tracking

#include "MiddleLayer.h"

#include <atlstr.h>

// Pass through evaluation which also performs some recording
float Middle_Layer::evaluate(const vector<bool>& solution) {
  // keeps track of how many times "evaluate" is called
  evaluations++;
  float fitness = evaluator->evaluate(solution);

  // keeps track of improvements in best found fitness
  if (best_fitness < fitness) {
	  
	float unitation = this->unitation(solution);
	time_t time = ::time(nullptr) - start_time;

	if (log)
	{
		CString log_message;
		log_message.AppendFormat("previous best fitness: %f", best_fitness);
		log->vPrintLine(log_message);
	}

    best_fitness = fitness;
    best_solution = solution;
    best_found = evaluations;
	best_unitation = unitation;

	if (log) {
		CString log_message;
		log_message.AppendFormat("best fitness: %f best unitation: %f ffe: %u time: %u", best_fitness, best_unitation, evaluations, (unsigned int)time);

		CString best;

		for (size_t i = 0; i < solution.size(); i++)
		{
			if (i % 10 == 0 && i > 0)
			{
				best.Append(" ");
			}

			best.AppendFormat("%d", solution.at(i) ? 1 : 0);
		}

		log->vPrintLine(log_message, true);
		log->vPrintLine(best);
	}

    results.add(fitness, unitation, evaluations, time, solution);
    // If configured to do so, will append the current best
    // to the output solution file
    if (!disable_solution_outfile) {
      std::ofstream output(outfile, std::fstream::app);
      // Format is <fitness> <evaluations> <solution>
      output << best_fitness << " " << evaluations << " ";
      print(best_solution, output);
    }
  }
  return fitness;
}

float Middle_Layer::unitation(const vector<bool>& solution) {
	int unitation = 0;

	for (size_t i = 0; i < solution.size(); i++) {
		if (solution.at(i) == 1) {
			unitation++;
		}
	}

	return (float)unitation / (float)solution.size();
}