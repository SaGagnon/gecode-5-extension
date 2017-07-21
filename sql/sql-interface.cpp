#include "sql-interface.hh"

#include <algorithm>
#include <sqlite3.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <unordered_map>

#include "sql-str-query.hpp"

//TODO: Faire une description ici qui mentionne que l'on doit toujours caller
//TODO: end_execution() sinon la DB finit dans un état incorrect + mentionner
//TODO: les attributs globales.

namespace CBSDB {

  /*****************************************************************************
   * BEGIN GLOBAL PROPERTIES
   ****************************************************************************/
  const int INVALID = -1;

  sqlite3 *current_db = nullptr;

  sqlite3_int64 current_exec_id = INVALID;
  int current_node_id = INVALID;

  EXEC_TYPE exec_type = EXEC_TYPE::INVALID_EXEC;

  bool unique_sol_found = false;

  std::unordered_map<unsigned int, int> correctVal;


  /*****************************************************************************
   * HELPER METHODS
   ****************************************************************************/

  bool db_active() {
    return current_db != nullptr;
  }

  bool insertion_exec_type() {
    assert(exec_type != EXEC_TYPE::INVALID_EXEC);
    return exec_type == EXEC_TYPE::INSERT_INTO_DB;
  }

  bool get_uniq_sol_exec_type() {
    assert(exec_type != EXEC_TYPE::INVALID_EXEC);
    return exec_type == EXEC_TYPE::GET_UNIQUE_SOLUTION;
  }

  void set_uniq_sol_found(bool b) {
    unique_sol_found = b;
  }

  void exit_on_error(const std::string& err = "");

  void db_exec(const std::string& sql_statement,
               int (*callback)(void*,int,char**,char**) = nullptr,
               void *first_arg_callback = nullptr) {
    char *errMsg;
    if (sqlite3_exec(current_db, sql_statement.c_str(),
                     callback, first_arg_callback, &errMsg) != SQLITE_OK) {
      std::cout << sql_statement << std::endl;
      std::cout << "err: " << errMsg << std::endl;
      sqlite3_free(errMsg);
      exit_on_error();
    }
  }

  void close_db() {
    if (!db_active()) return;
    db_exec("END TRANSACTION;");
    sqlite3_close(current_db);
    current_db = nullptr;
  }

  void exit_on_error(const std::string& err) {
    std::cout << err << std::endl;
    close_db();
    exit(-1);
  }

  int get_exec_id(void* exec_id, int argc, char** argv, char** colName) {
    if (argv[0] == nullptr)
      *(sqlite3_int64*)exec_id = INVALID;

    *(sqlite3_int64*)exec_id = std::stoi(argv[0]);
    return 0;
  }

  int get_exec_results(void* notUsed, int argc, char** argv, char** colName) {
    auto var_id = (unsigned int) std::stoi(argv[0]);
    int val = std::stoi(argv[1]);
    correctVal[var_id] = val;
    return 0;
  }

  /*****************************************************************************
   * IMPLEMENTATION
   ****************************************************************************/

  void set_exec_type(EXEC_TYPE e) {
    exec_type = e;
  }

  void start_execution(std::string pb_name, unsigned int num_ex,
                       std::string branching_name, std::string path_db) {

    if (exec_type == EXEC_TYPE::INVALID_EXEC)
      exit_on_error("Execution type is invalid");

    if (sqlite3_open(path_db.c_str(), &current_db) != SQLITE_OK)
      exit_on_error("Error in opening database");

    db_exec("BEGIN TRANSACTION;");

    if (exec_type == EXEC_TYPE::GET_UNIQUE_SOLUTION) {
      db_exec(sql_str_insert_into_executions(
          {0, 0, pb_name, num_ex, branching_name}
      ));
      current_exec_id = sqlite3_last_insert_rowid(current_db);
      assert(current_exec_id != INVALID);
    } else if (exec_type == EXEC_TYPE::INSERT_INTO_DB) {
      {
        std::stringstream sql;
        sql <<
            "SELECT exec_id FROM executions "
              "WHERE pb_name = '" << pb_name << "' "
              "AND num_ex = " << num_ex << " "
              "AND single_sol_found = 1 "
              "AND nodes_inserted = 0;";
        db_exec(sql.str(), get_exec_id, &current_exec_id);
      }

      if (current_exec_id == INVALID)
        exit_on_error("No single solution found for given problem and num_ex");

      {
        std::stringstream sql;
        sql << "SELECT var_id, val FROM results "
          "WHERE exec_id = " << current_exec_id << ";";
        db_exec(sql.str(), get_exec_results);
      }
    }
  }

  void new_node() {
    if (!db_active()) return;
    if (!insertion_exec_type()) return;
    current_node_id++;

    db_exec(sql_str_insert_into_nodes(
      {(unsigned int) current_exec_id, (unsigned int)current_node_id}
    ));
  }

  void new_propagator(unsigned int prop_id) {
    if (!db_active()) return;
    if (!insertion_exec_type()) return;

    db_exec(sql_str_insert_into_propagators(
      {(unsigned int)current_exec_id, (unsigned int)current_node_id, prop_id}
    ));
  }


  void insert_varval_density(unsigned int prop_id, unsigned int var_id, int val,
                             double max_sd, double a_avg_sd,
                             unsigned int var_dom_size, double max_rel_sd,
                             double max_rel_ratio) {
    if (!db_active()) return;
    if (!insertion_exec_type()) return;

    db_exec(sql_str_insert_into_densities({
          (unsigned int)current_exec_id,
          (unsigned int)current_node_id,
          prop_id,
          var_id,
          val,
          max_sd,
          a_avg_sd,
          var_dom_size,
          max_rel_sd,
          max_rel_ratio
    }));
  }

  void insert_varval_in_solution(unsigned int var_id, int val) {
    if (!db_active()) return;
    if (!get_uniq_sol_exec_type()) return;

    db_exec(sql_str_insert_into_results(
      {(unsigned int)current_exec_id, var_id, val}
    ));
  }

  void end_execution() {
    if (!db_active()) return;

    auto update_executions_set_1 = [&](std::string param) {
      std::stringstream sql;
      sql <<"UPDATE executions SET " << param << " = 1 "
            "WHERE exec_id = " << current_exec_id << ";";
      db_exec(sql.str());
    };

    if (exec_type == EXEC_TYPE::GET_UNIQUE_SOLUTION && unique_sol_found)
      update_executions_set_1("single_sol_found");
    else if (exec_type == EXEC_TYPE::INSERT_INTO_DB)
      update_executions_set_1("nodes_inserted");

    close_db();
  }
}

