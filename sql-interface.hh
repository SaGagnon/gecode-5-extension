#ifndef __SQL_INTERFACE_HH__
#define __SQL_INTERFACE_HH__

#include <string>
#include <utility>
#include <vector>

#define CBSDB_FAILED -1
#define CBSDB_SUCCESS  0
#define CBSDB_NO_ACTION_TAKEN 1


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
   *
   * @param problem_name
   *  The name of the problem we are solving.
   * @param num_example
   *  An id that differentiate different examples of the same problem.
   * @param solver_name
   *  Name of the CP Solver.
   * @param ech_method
   *  Method that must be use for sampling the nodes during the execution.
   * @param max_nb_nodes
   *  Maximum number of nodes to insert into the database.
   * @param branching_heuristic_name
   *  Name of the branching heuristic used in the model for the problem.
   * @param model_version_name
   *  This parameter can be used to identify different version of CP models for
   *  the problem.
   * @param path_to_db
   *  Path to file containing the database we want to use.
   * @return
   *  CBSDB_SUCCESS: Database successfully opened and execution label created.
   *  CBSDB_FAILED: Error in opening database or creating execution label.
   */
  int start_execution(std::string problem_name, unsigned int num_example,
                      std::string solver_name, EchMethod ech_method,
                      unsigned int max_nb_nodes,
                      std::string branching_heuristic_name,
                      std::string model_version_name, std::string path_to_db);

  /**
   * Insert a new node associated with the current execution. By default, every
   * new node is considered to be satisfiable node (all the variable assigned
   * are part of the solution). We update this information in the method
   * "end_execution".
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @return
   *  CBSDB_SUCCESS: Node successfully inserted in database.
   *  CBSDB_FAILED: Failed to insert node in database.
   *  CBSDB_NO_ACTION_TAKEN: We didn't try to insert the node According the
   *  sampling method given in start_execution, we must skip this one.
   */
  int new_node();

  /**
   * Insert a new propagator associated with the current node.
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @param propagator_name
   *  Name of the propagator.
   * @param consistency_level_name
   *  Name of the consistency level of the propagator
   * @param solution_count
   *  Number of possible solutions in the propagator.
   * @return
   *  CBSDB_SUCCESS: Propagator successfully inserted in database.
   *  CBSDB_FAILED: Failed to insert propagator in database.
   *  CBSDB_NO_ACTION_TAKEN: We didn't try to insert the propagator. According
   *  the sampling method given in start_execution, we must skip this one.
   */
  int new_propagator(std::string propagator_name, unsigned int prop_id,
                     std::string consistency_level_name,
                     double solution_count);

  /**
   * Insert a new (var,val) pair associated with the current propagator.
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @param var_idx
   *  An ID that uniquly defines the corresponding variable during the whole
   *  execution.
   * @param val
   *  A value in the domain of the variable.
   * @param dens
   *  The calculated density associated with the pair (var,val).
   */
  int insert_varval_density(unsigned int prop_id, unsigned int var_id,
                            int val, double dens);


  int insert_varval_density_features(
    unsigned int prop_id, unsigned int var_id, int val, double dens,
    double sln_cnt, double sum_sln_cnt, double a_avg_sd, double var_dom_size,
    double var_dens_entropy, double max_rel_sd, double max_rel_ratio,
    double w_sc_avg, double w_anti_sc_avg);

  /**
   * Create a new solution linked to the current execution.
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @return
   *  CBSDB_SUCCESS: New solution successfully created.
   *  CBSDB_NO_ACTION_TAKEN: See method description.
   */
  int new_solution();

  /**
   *
   * Insert a pair in the last solution created.
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @param var_idx
   *  An ID that uniquly defines the corresponding variable during the whole
   *  execution.
   * @param val
   *  A value in the domain of the variable that is part of the current solution.
   * @return
   *  CBSDB_SUCCESS: (var,val) pair successfully inserted.
   *  CBSDB_FAILED: Error inserting (var,val) pair.
   *  CBSDB_NO_ACTION_TAKEN: See method description.
   */
  int insert_varval_in_solution(unsigned int var_idx, int val);

  /**
   * Close the database. We also update all nodes that contain a pair (var,val)
   * with a density of 100% which is not part of the solution to be
   * unsatifiable.
   *
   * If start_execution wasn't called, the methond will do nothing and return
   * CBSD_NO_ACTION_TAKEN.
   *
   * @return
   *  CBSDB_SUCCESS: Database successfully closed.
   *  CBSDB_FAILED: Error when closing database.
   *  CBSDB_NO_ACTION_TAKEN: See method description.
   */
  int end_execution();
}


/**
 * Gecode extension
 */
#include <gecode/int.hh>
namespace CBSDB {
  int insert_if_solution(const Gecode::IntVarArray& x);
}

#endif
