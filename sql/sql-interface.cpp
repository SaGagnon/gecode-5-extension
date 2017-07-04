#include "sql-interface.hh"

#include <algorithm>
#include <sqlite3.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>

#include "sql-str-query.hpp"

//TODO: Faire une description ici qui mentionne que l'on doit toujours caller
//TODO: end_execution() sinon la DB finit dans un état incorrect + mentionner
//TODO: les attributs globales.

namespace CBSDB {

  /*****************************************************************************
   * BEGIN GLOBAL PROPERTIES
   ****************************************************************************/
  const int INVALID = -1;

  sqlite3 *current_db = NULL;
  unsigned int ech_method;

  sqlite3_int64 current_exec_id = INVALID;
  int current_node_id = INVALID;
  int current_result_id = INVALID;

  bool solution_found = false;

  unsigned int max_nb_nodes;

  /*****************************************************************************
   * HELPER METHODS
   ****************************************************************************/
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


  bool db_exec(std::string sql_statement) {
    char *errMsg;
    if (sqlite3_exec(current_db, sql_statement.c_str(),
                     NULL, NULL, &errMsg) != SQLITE_OK) {
      std::cout << sql_statement;
      std::cout << "err: " << errMsg << std::endl;
      sqlite3_free(errMsg);
      return false;
    }
    return true;
  }

  /*****************************************************************************
   * DEFINES
   ****************************************************************************/

#define EXEC_SQL(str) \
  if (!db_exec(str)) \
    return;

#define DB_ACTIVE_OR_RETURN if (current_db == NULL) return;

#define CURRENT_EXECID_VALID_OR_RETURN \
  if (current_exec_id == INVALID) { \
    std::cout << "Current execution ID invalid" << std::endl; \
    return; \
  }

#define DO_WE_INSERT_OR_RETURN if(!do_we_insert()) return;

  /*****************************************************************************
   * IMPLEMENTATION
   ****************************************************************************/

  void start_execution(struct executions& s, std::string path_to_db) {
    if (sqlite3_open(path_to_db.c_str(), &current_db) != SQLITE_OK) {
      std::cout << "Error in opening database" << std::endl;
      return;
    }

    if (s.ech_method < 0 || s.ech_method >= NB_ECH_METHOD) {
      std::cout << "Sampling method not valid." << std::endl;
      return;
    }

    EXEC_SQL("BEGIN TRANSACTION;")

    ech_method = s.ech_method;
    max_nb_nodes = s.max_nb_nodes;
    EXEC_SQL(sql_str_insert_into_executions(s))

    // ID of the executions we just inserted. Will be used for inserting nodes.
    current_exec_id = sqlite3_last_insert_rowid(current_db);
  }

  void new_node() {
    DB_ACTIVE_OR_RETURN
    current_node_id++;
    DO_WE_INSERT_OR_RETURN
    CURRENT_EXECID_VALID_OR_RETURN

    nodes n;
    n.exec_id = (unsigned int)current_exec_id;
    n.node_id = (unsigned int)current_node_id;
    n.sat = 0;

    EXEC_SQL(sql_str_insert_into_nodes(n))
  }

  void insert_varval_density(struct densities& s) {
    DB_ACTIVE_OR_RETURN
    DO_WE_INSERT_OR_RETURN

    EXEC_SQL(sql_str_insert_into_densities(s))
  }

  void insert_varval_in_assigned(struct assigned& s) {
    DB_ACTIVE_OR_RETURN
    DO_WE_INSERT_OR_RETURN

    EXEC_SQL(sql_str_insert_into_assigned(s))
  }

  void new_solution() {
    DB_ACTIVE_OR_RETURN

    current_result_id++;
  }

  void insert_varval_in_solution(struct results& s) {
    DB_ACTIVE_OR_RETURN
    CURRENT_EXECID_VALID_OR_RETURN

    EXEC_SQL(sql_str_insert_into_results(s))
    solution_found = true;
  }

  void end_execution() {
    DB_ACTIVE_OR_RETURN

    // We find unsatifiable nodes and update them in the database.
    if(solution_found) {
      std::stringstream sql;
      sql << " UPDATE nodes SET sat = 1"
          << " WHERE exec_id = " << current_exec_id
          << " AND node_id IN ("

        << " SELECT nn.node_id"
        << " FROM nodes AS nn"
        << " WHERE"
        << "       ( SELECT count(*)"
        << "         FROM assigned AS a"
        << "           JOIN results AS r"
        << "             ON a.exec_id = r.exec_id"
        << "                AND a.var_idx = r.var_idx"
        << "                AND a.val = r.val"
        << "                AND r.res_id = 0" //TODO: temporaire
        << "         WHERE a.exec_id = nn.exec_id"
        << "               AND a.node_id = nn.node_id"
        << "       ) == ("
        << "         SELECT count(*)"
        << "         FROM assigned AS a"
        << "         WHERE a.exec_id = nn.exec_id"
        << "               AND a.node_id = nn.node_id"
        << "       )"
        << " AND nn.exec_id = " << current_exec_id

        << ");";

      EXEC_SQL(sql.str())
    }

    EXEC_SQL("END TRANSACTION;")
    sqlite3_close(current_db);
    current_db = NULL;
  }

}

namespace CBSDB {
//  void insert_if_solution(const Gecode::IntVarArray& x) {
//    DB_ACTIVE_OR_RETURN
//    for (int i=0; i<x.size(); i++)
//      if (!x[i].assigned())
//        return;
//
//    CBSDB::new_solution();
//    for(int i=0; i<x.size(); i++) {
//      CBSDB::insert_varval_in_solution(x[i].varimp()->id(), x[i].val());
//    }
//  }

//  void insert_if_solution(const Gecode::BoolVarArray& x) {
//    DB_ACTIVE_OR_RETURN
//    for (int i=0; i<x.size(); i++)
//      if (!x[i].assigned())
//        return;
//
//    int ret = CBSDB::new_solution();
//    if (ret) return ret;
//    for(int i=0; i<x.size(); i++) {
//      ret = CBSDB::insert_varval_in_solution(x[i].varimp()->id(), x[i].val());
//      if (ret) return ret;
//    }
//  }
}
