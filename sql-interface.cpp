#include "sql-interface.hh"

#include <algorithm>
#include <sqlite3.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>

//TODO: Faire une description ici qui mentionne que l'on doit toujours caller
//TODO: end_execution() sinon la DB finit dans un état incorrect + mentionner
//TODO: les attributs globales.

namespace CBSDB {

  const int INVALID = -1;

  // GLOBAL PROPERTIES
  sqlite3 *current_db = NULL;
  EchMethod ech_method;

  sqlite3_int64 current_exec_id = INVALID;
  int current_node_id = INVALID;
  int current_result_id = INVALID;

  bool solution_found = false;

  unsigned int max_nb_nodes;


  std::string to_lower_case(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower;
  }


  bool do_we_insert() {
    if (current_node_id >= max_nb_nodes)
      return false;

    switch (ech_method) {
      case FULL:
        return true;
      case ECH10:
        return current_node_id % 10 == 0;
      case ECH100:
        return current_node_id % 100 == 0;
      case ECH1000:
        return current_node_id % 1000 == 0;
      default:
        assert(false);
    }
    assert(false);
  }


  int db_exec(std::string sql_statement, std::string err) {
    char *errMsg;
    if (sqlite3_exec(current_db, sql_statement.c_str(),
                     NULL, NULL, &errMsg) != SQLITE_OK) {
      std::cout << err << ": " << errMsg << std::endl;
      sqlite3_free(errMsg);
      return CBSDB_FAILED;
    }
    return CBSDB_SUCCESS;
  }


  int start_execution(std::string problem_name, unsigned int num_example,
                      std::string solver_name, EchMethod ech_method0,
                      unsigned int max_nb_nodes0,
                      std::string branching_heuristic_name,
                      std::string model_version_name, std::string path_to_db) {
    if (sqlite3_open(path_to_db.c_str(), &current_db) != SQLITE_OK) {
      std::cout << "Error in opening database" << std::endl;
      return CBSDB_FAILED;
    }

    if (ech_method0 < 0 || ech_method0 >= NB_ECH_METHOD) {
      std::cout << "Sampling method not valid." << std::endl;
      return CBSDB_FAILED;
    }
    ech_method = ech_method0;
    max_nb_nodes = max_nb_nodes0;

    if (db_exec("BEGIN TRANSACTION;", "Begin transaction failed"))
      return CBSDB_FAILED;

    std::stringstream sql;
    sql
      << "insert into executions("
      << "pb_name, num_ex, solveur_name, ech_method, max_nb_nodes, "
      <<  "branching_name, version_name) "
      << "values("
      << "'" << to_lower_case(problem_name) << "', "
      << num_example << ", "
      << "'" << to_lower_case(solver_name) << "', "
      << ech_method0 << ", "
      << max_nb_nodes0 << ", "
      << "'" << to_lower_case(branching_heuristic_name) << "', "
      << "'" << to_lower_case(model_version_name) << "'); ";

    if (db_exec(sql.str(), "Creation of new execution failed"))
      return CBSDB_FAILED;

    // ID of the executions we just inserted. Will be used for inserting nodes.
    current_exec_id = sqlite3_last_insert_rowid(current_db);
    return CBSDB_SUCCESS;
  }


  int new_node() {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;
    current_node_id++;

    if (do_we_insert()) {
      if (current_exec_id == INVALID) {
        std::cout << "Current execution ID invalid" << std::endl;
        return CBSDB_FAILED;
      }

      std::stringstream sql;
      sql << "insert into nodes(node_id, exec_id, sat) values("
          << current_node_id << ", "
          << current_exec_id << ","
          << 0 << ");";

      if (db_exec(sql.str(), "Creation of new node failed"))
        return CBSDB_FAILED;

      return CBSDB_SUCCESS;
    }
    return CBSDB_NO_ACTION_TAKEN;
  }


  int new_propagator(std::string propagator_name, unsigned int prop_id,
                     std::string consistency_level_name,
                     double solution_count) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;

    if (do_we_insert()) {
      if (current_node_id == INVALID) {
        std::cout << "Current node ID invalid" << std::endl;
        return CBSDB_FAILED;
      }

      //TODO: Hack avec propagator_name.
      std::stringstream sql;
      sql << "insert into propagators(exec_id, node_id, prop_id, prop_name, "
          << "cons_lvl, log_solutionCount) values("
          << current_exec_id << ", "
          << current_node_id << ", "
          << prop_id << ", "
          << "'" << propagator_name << "', "
          << "'" << consistency_level_name << "', "
          << log(solution_count) << ");";

      if (db_exec(sql.str(), "Creation of new propagtor failed"))
        return CBSDB_FAILED;

      return CBSDB_SUCCESS;
    }
    return CBSDB_NO_ACTION_TAKEN;
  }

  int insert_varval_density(unsigned int prop_id, unsigned int var_id, int val,
                            double dens) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;

    if (do_we_insert()) {
      std::stringstream sql;
      sql << "insert into densities(exec_id, node_id, prop_id, var_idx, val, "
        "dens) values(";
      sql << current_exec_id << ","
          << current_node_id << ","
          << prop_id << ","
          << var_id << ","
          << val << ","
          // We encode the density value between 0 and 240 to save space.
          << (int) (dens * 240) << ");";

      if (db_exec(sql.str(), "Creation of the (var,val) pair failed"))
        return CBSDB_FAILED;
    }
    return CBSDB_NO_ACTION_TAKEN;
  }


  int insert_varval_density_features(
    unsigned int prop_id, unsigned int var_id, int val, double dens,
    double sln_cnt, double sum_sln_cnt, double a_avg_sd, double var_dom_size,
    double var_dens_entropy, double max_rel_sd, double max_rel_ratio,
    double w_sc_avg, double w_anti_sc_avg, double w_t_avg, double w_anti_t_avg,
    double w_d_avg) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;

    if (do_we_insert()) {
      std::stringstream sql;
      sql << "insert into densities(exec_id, node_id, prop_id, var_idx, "
        "val, dens, log_sln_cnt, log_sum_sln_cnt, a_avg_sd, var_dom_size, "
        "var_dens_entropy, max_rel_sd, max_rel_ratio, w_sc_avg, "
        "w_anti_sc_avg, w_t_avg, w_anti_t_avg, w_d_avg) values(";
      sql << current_exec_id << ","
          << current_node_id << ","
          << prop_id << ","
          << var_id << ","
          << val << ","
          << dens << ","
          << log(sln_cnt) << ","
          << log(sum_sln_cnt) << ","
          << a_avg_sd << ","
          << var_dom_size << ","
          << var_dens_entropy << ","
          << max_rel_sd << ","
          << max_rel_ratio << ","
          << w_sc_avg << ","
          << w_anti_sc_avg << ","
          << w_t_avg << ","
          << w_anti_t_avg << ","
          << w_d_avg << ");";

      if (db_exec(sql.str(), "Creation of the (var,val) pair failed"))
        return CBSDB_FAILED;
    }
    return CBSDB_NO_ACTION_TAKEN;

  }


  int new_solution() {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;
    current_result_id++;
    return CBSDB_SUCCESS;
  }

  int insert_varval_in_solution(unsigned int var_idx, int val) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;

    if (current_exec_id == INVALID) {
      std::cout << "Current execution ID invalid" << std::endl;
      return CBSDB_FAILED;
    }

    std::stringstream sql;
    sql << "insert into results(res_id, exec_id, var_idx, val) "
        << "values("
        << current_result_id << ", "
        << current_exec_id << ", "
        << var_idx << ", "
        << val << ");";

    if (db_exec(sql.str(), "Insertion in results failed"))
      return CBSDB_FAILED;

    solution_found = true;
    return CBSDB_SUCCESS;
  }


  int end_execution() {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;

    // We find unsatifiable nodes and update them in the database.
    if(solution_found) {
      std::stringstream sql;
      sql << "UPDATE nodes SET sat = 1"
          << " WHERE exec_id = " << current_exec_id
          << " AND node_id IN ("

        << " SELECT nn.node_id"
        << " FROM nodes AS nn"
        << " WHERE"
        << " NOT exists("
        << "     SELECT *"
        << "     FROM nodes AS n"
        << "       JOIN densities AS d"
        << "         ON n.exec_id = d.exec_id"
        << "            AND n.node_id = d.node_id"
        << "       LEFT JOIN results AS r"
        << "         ON d.exec_id = r.exec_id"
        << "            AND d.var_idx = r.var_idx"
        << "            AND d.val = r.val"
        << "            AND r.res_id = 0" //TODO: temporaire
        << "     WHERE n.exec_id = nn.exec_id AND n.node_id = nn.node_id"
        << "           AND d.dens = 1"
        << "           AND r.exec_id IS NULL"
        << " ) AND ("
        << "         SELECT count(DISTINCT d.var_idx)"
        << "         FROM nodes AS n"
        << "           JOIN densities AS d"
        << "             ON n.exec_id = d.exec_id"
        << "                AND n.node_id = d.node_id"
        << "           JOIN results AS r"
        << "             ON d.exec_id = r.exec_id"
        << "                AND d.var_idx = r.var_idx"
        << "                AND d.val = r.val"
        << "                AND r.res_id = 0" //TODO: temporaire
        << "         WHERE n.exec_id = nn.exec_id"
        << "               AND n.node_id = nn.node_id"
        << "       ) == ("
        << "         SELECT count(DISTINCT d.var_idx)"
        << "         FROM nodes AS n"
        << "           JOIN densities AS d"
        << "             ON n.exec_id = d.exec_id"
        << "                AND n.node_id = d.node_id"
        << "         WHERE n.exec_id = nn.exec_id"
        << "               AND n.node_id = nn.node_id"
        << "       )"
        << " AND nn.exec_id = " << current_exec_id

        << ");";


      if (db_exec(sql.str(), "Update of unsatisfiable nodes failed"))
        return CBSDB_FAILED;
    }

    if (db_exec("END TRANSACTION;", "End transaction failed"))
      return CBSDB_FAILED;
    sqlite3_close(current_db);
    current_db = NULL;
    return CBSDB_SUCCESS;
  }

}

namespace CBSDB {
  int insert_if_solution(const Gecode::IntVarArray& x) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;
    for (int i=0; i<x.size(); i++)
      if (!x[i].assigned())
        return CBSDB_NO_ACTION_TAKEN;

    int ret = CBSDB::new_solution();
    if (ret) return ret;
    for(int i=0; i<x.size(); i++) {
      ret = CBSDB::insert_varval_in_solution(x[i].varimp()->id(), x[i].val());
      if (ret) return ret;
    }
    return CBSDB_SUCCESS;
  }

  int insert_if_solution(const Gecode::BoolVarArray& x) {
    if (current_db == NULL) return CBSDB_NO_ACTION_TAKEN;
    for (int i=0; i<x.size(); i++)
      if (!x[i].assigned())
        return CBSDB_NO_ACTION_TAKEN;

    int ret = CBSDB::new_solution();
    if (ret) return ret;
    for(int i=0; i<x.size(); i++) {
      ret = CBSDB::insert_varval_in_solution(x[i].varimp()->id(), x[i].val());
      if (ret) return ret;
    }
    return CBSDB_SUCCESS;
  }
}
