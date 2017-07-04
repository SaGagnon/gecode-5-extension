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

  enum EchMethod {
    FULL,
    ECH10,
    ECH100,
    ECH1000,
    NB_ECH_METHOD
  };

  /**
   * This function must be called first before the solver begins to solve its
   * model. Its role is to connect to the database and insert a new execution
   * label in the database. The database must have been previously initialized
   * with the correct tables.
   */
  void start_execution(struct executions&, std::string path_to_db);

  /**
   * Insert a new node associated with the current execution. By default, every
   * new node is considered to be satisfiable node (all the variable assigned
   * are part of the solution). We update this information in the method
   * "end_execution".
   */
  void new_node();

  /**
   * Insert a new (var,val) pair associated with the current node in densities
   * table.
   */
  void insert_varval_density(struct densities&);


  /**
   * Insert a new (var,val) pair associated with the current node in assigned
   * table.
   */
  void insert_varval_in_assigned(struct assigned&);

  /**
   * Create a new solution linked to the current execution.
   */
  void new_solution();

  /**
   * Insert a pair in the last solution created.
   */
  void insert_varval_in_solution(struct results&);

  /**
   * Close the database. We also update all nodes that contain a pair (var,val)
   * with a density of 100% which is not part of the solution to be
   * unsatifiable.
   */
  void end_execution();
}

/**
 * Gecode extension
 */
#include <gecode/int.hh>
namespace CBSDB {
  void insert_if_solution(const Gecode::IntVarArray& x);
  void insert_if_solution(const Gecode::BoolVarArray& x);
}

#endif
