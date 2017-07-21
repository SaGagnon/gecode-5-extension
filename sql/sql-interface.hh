#ifndef __SQL_INTERFACE_HH__
#define __SQL_INTERFACE_HH__

#include <string>
#include <utility>
#include <vector>

#include "sql-structs.hh"

/**
 * Generic interface
 */
namespace CBSDB {

  enum EXEC_TYPE {
    GET_UNIQUE_SOLUTION,
    INSERT_INTO_DB,
    LEN,
    INVALID_EXEC
  };

  // TODO: Commentaire...
  void set_exec_type(EXEC_TYPE);

  /**
   * This function must be called first before the solver begins to solve its
   * model. Its role is to connect to the database and insert a new execution
   * label in the database. The database must have been previously initialized
   * with the correct tables.
   */
  void start_execution(std::string pb_name, unsigned int num_ex,
                       std::string branching_name, std::string path_db);

  /**
   * Insert a new node associated with the current execution.
   */
  void new_node();

  /**
   * Insert a new propagator associated with the current node.
   */
  void new_propagator(unsigned int prop_id);

  /**
   * Insert a new (var,val) pair associated with the current node.
   */
  void insert_varval_density(unsigned int prop_id, unsigned int var_id, int val,
                             double max_sd, double a_avg_sd,
                             unsigned int var_dom_size, double max_rel_sd,
                             double max_rel_ratio);

  /**
   * Insert correct assigment in solution.
   */
  void insert_varval_in_solution(unsigned int var_id, int val);

  void end_execution();
}

/**
 * Gecode extension
 */
#include <gecode/int.hh>
#include <unordered_map>

namespace CBSDB {

  bool db_active();
  bool insertion_exec_type();
  bool get_uniq_sol_exec_type();
  void set_uniq_sol_found(bool);

  extern std::unordered_map<unsigned int, int> correctVal;

  template<class Arr>
  bool is_node_sat(const Arr& x) {
    for (int i=0; i<x.size(); i++)
      if (x[i].assigned() && x[i].val() != correctVal[x[i].id()])
        return false;
    return true;
  }

  template<class Arr>
  void insert_if_solution(const Arr& x) {
    if (!db_active()) return;
    if (!get_uniq_sol_exec_type()) return;

    for (int i=0; i<x.size(); i++)
      if (!x[i].assigned())
        return;

    set_uniq_sol_found(true);

    for(int i=0; i<x.size(); i++)
      insert_varval_in_solution(x[i].varimp()->id(), x[i].val());
  }
}

#endif
