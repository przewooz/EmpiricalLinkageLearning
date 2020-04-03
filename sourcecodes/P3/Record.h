// Brian Goldman

// A record stores the result information for a single optimization run,
// recording how many evaluations are required to reach each new fitness level

#ifndef RECORD_H_
#define RECORD_H_

#include <vector>
#include <utility>
#include <array>
#include <stdlib.h>
#include "Util.h"
#include "Configuration.h"

using std::vector;
using std::tuple;
using std::array;

#include <ctime>

// Labels used for decoding the "summary" return statistics
enum RecordKey {
  MES,
  MAD,
  FAILURES
};

// Stores optimization history information
class Record {
 public:
  Record(Configuration& config)
      : verbosity(config.get<int>("verbosity")) {
  }
  virtual ~Record() = default;
  // Tracks increases in fitness.  If verbosity is set high enough, outputs
  // information to the screen.
  void add(float fitness, float unitation, int evals, time_t time, const vector<bool>& solution);
  // The best fitness achieved in this optimization run
  const tuple<float, float, int, time_t, vector<bool>>& best() const;
  // Summarizes the end of run statistics for mutliple records
  static array<int, 3> summarize(const vector<Record>& records, Configuration& config);
  // access to the stored values:
  const vector<tuple<float, float, int, time_t, vector<bool>>>& progression() {
    return levels;
  }

 private:
  // Raw stored data
  vector<tuple<float, float, int, time_t, vector<bool>>> levels;

  int verbosity;

};

#endif /* RECORD_H_ */
