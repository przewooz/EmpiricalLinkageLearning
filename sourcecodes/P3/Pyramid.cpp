/*
 * Pyramid.cpp
 *
 *  Created on: Jun 26, 2013
 *      Author: goldman
 */

#include "Pyramid.h"

// Applies crossover between the passed in solution as each level
// of the pyramid
bool Pyramid::climb(vector<bool> & solution, float & fitness) {
	return climb(solution, fitness, 0, false, false);
}

bool Pyramid::climb(vector<bool> & solution, float & fitness, vector<size_t> & addedToLevels) {
	bool b_added = false;

	if (add_unique(solution, 0)) {
		addedToLevels.push_back(0);
		b_added = true;
	}

	for (size_t lev = 0; lev < pops.size(); lev++) {
		float prev = fitness;
		// Use population clusters and population solutions to make new solution
		pops[lev].improve(rand, solution, fitness, evaluator);
		// add it to the next level if its a strict fitness improvement,
		// or configured to always add solutions
		if (!only_add_improvements || prev < fitness) {
			if (add_unique(solution, lev + 1)) {
				addedToLevels.push_back(lev + 1);
				b_added = true;
			}
		}
	}

	return b_added;
}

bool Pyramid::climb(vector<bool> & solution, float & fitness, size_t level, bool hillClimber, bool climbOnlyIfUnique) {
	if (hillClimber) {
		hill_climber(rand, solution, fitness, evaluator);
	}

	// attempts to add this solution to the base level
	bool b_added = add_unique(solution, level);

	if (!climbOnlyIfUnique || b_added) {
		for (size_t lev = level; lev < pops.size(); lev++) {
			float prev = fitness;
			// Use population clusters and population solutions to make new solution
			pops[lev].improve(rand, solution, fitness, evaluator);
			// add it to the next level if its a strict fitness improvement,
			// or configured to always add solutions
			if (!only_add_improvements || prev < fitness) {
				b_added = add_unique(solution, lev + 1) || b_added;
			}
		}
	}

	return b_added;
}

bool Pyramid::feedback(vector<bool> & solution, float & fitness) {
	bool b_added = false;
	
	for (size_t lev = 0; lev < pops.size(); lev++) {
		float prev = fitness;
		// Use population clusters and population solutions to make new solution
		pops[lev].improve(rand, solution, fitness, evaluator);
		// add it to the next level if its a strict fitness improvement,
		// or configured to always add solutions
		if (prev < fitness) {
			b_added = add_unique(solution, pops.size()) || b_added;
		}
	}

	return b_added;
}

// Add the solution to the specified level of the population, only if it is
// unique.  Returns true if added.
bool Pyramid::add_unique(const vector<bool> & solution, size_t level) {
  if (seen.find(solution) == seen.end()) {
    if (pops.size() == level) {
      // Create new levels as necessary
      pops.push_back(Population(config, evaluator));
    }
    // Add the solution and rebuild the tree
    pops[level].add(solution);
    pops[level].rebuild_tree(rand);
    // track that this solution is now in the pyramid
    seen.insert(solution);
    return true;
  }
  return false;
}

// Performs a full iteration of P3
bool Pyramid::iterate() {
  // generate a random solution
  vector<bool> solution = rand_vector(rand, length);
  float fitness = evaluator.evaluate(solution);
  // perform a local search hill climber
  hill_climber(rand, solution, fitness, evaluator);
  // perform crossover with each level of the pyramid
  climb(solution, fitness);
  // P3 never "converges"

  return true;
}
