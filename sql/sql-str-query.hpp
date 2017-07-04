// This file is autogenerated

namespace CBSDB {

std::string sql_str_insert_into_executions(struct executions& s) {
std::stringstream sql;
sql
<< "insert into executions ("
<< "exec_id," 
<< "pb_name," 
<< "num_ex," 
<< "ech_method," 
<< "max_nb_nodes," 
<< "branching_name," 
<< ") values ("
<< s.exec_id << ", "
<< "'" << s.pb_name << "',"
<< s.num_ex << ", "
<< s.ech_method << ", "
<< s.max_nb_nodes << ", "
<< "'" << s.branching_name << "',"
<< ");";
return str.str();
}

std::string sql_str_insert_into_results(struct results& s) {
std::stringstream sql;
sql
<< "insert into results ("
<< "exec_id," 
<< "res_id," 
<< "var_id," 
<< "val," 
<< ") values ("
<< s.exec_id << ", "
<< s.res_id << ", "
<< s.var_id << ", "
<< s.val << ", "
<< ");";
return str.str();
}

std::string sql_str_insert_into_nodes(struct nodes& s) {
std::stringstream sql;
sql
<< "insert into nodes ("
<< "exec_id," 
<< "node_id," 
<< "sat," 
<< ") values ("
<< s.exec_id << ", "
<< s.node_id << ", "
<< s.sat << ", "
<< ");";
return str.str();
}

std::string sql_str_insert_into_assigned(struct assigned& s) {
std::stringstream sql;
sql
<< "insert into assigned ("
<< "exec_id," 
<< "node_id," 
<< "var_id," 
<< "val," 
<< ") values ("
<< s.exec_id << ", "
<< s.node_id << ", "
<< s.var_id << ", "
<< s.val << ", "
<< ");";
return str.str();
}

std::string sql_str_insert_into_densities(struct densities& s) {
std::stringstream sql;
sql
<< "insert into densities ("
<< "exec_id," 
<< "node_id," 
<< "var_id," 
<< "val," 
<< "dens," 
<< "a_avg_sd," 
<< "var_dom_size," 
<< "max_rel_sd," 
<< "max_rel_ratio," 
<< "w_sc_avg," 
<< "w_anti_sc_avg," 
<< "w_t_avg," 
<< "w_anti_t_avg," 
<< "w_d_avg," 
<< ") values ("
<< s.exec_id << ", "
<< s.node_id << ", "
<< s.var_id << ", "
<< s.val << ", "
<< s.dens << ", "
<< s.a_avg_sd << ", "
<< s.var_dom_size << ", "
<< s.max_rel_sd << ", "
<< s.max_rel_ratio << ", "
<< s.w_sc_avg << ", "
<< s.w_anti_sc_avg << ", "
<< s.w_t_avg << ", "
<< s.w_anti_t_avg << ", "
<< s.w_d_avg << ", "
<< ");";
return str.str();
}

}
