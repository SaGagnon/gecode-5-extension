// This file is autogenerated

#ifndef __SQL_STRUCTS_HH__
#define __SQL_STRUCTS_HH__

namespace CBSDB {

 struct executions {
   unsigned int exec_id;
   std::string pb_name;
   unsigned int num_ex;
   unsigned int ech_method;
   unsigned int max_nb_nodes;
   std::string branching_name;
 };
 
 struct results {
   unsigned int exec_id;
   unsigned int res_id;
   unsigned int var_id;
   int val;
 };
 
 
 struct nodes {
   unsigned int exec_id;
   unsigned int node_id;
   unsigned int sat;
 };
 
 
 struct assigned {
   unsigned int exec_id;
   unsigned int node_id;
   unsigned int var_id;
   int val;
 };
 
 struct densities {
   unsigned int exec_id;
   unsigned int node_id;
   unsigned int var_id;
   int val;
   double dens;
   double a_avg_sd;
   unsigned int var_dom_size;
   double max_rel_sd;
   double max_rel_ratio;
   double w_sc_avg;
   double w_anti_sc_avg;
   double w_t_avg;
   double w_anti_t_avg;
   double w_d_avg;
 };

}

#endif //end __SQL_STRUCTS_HH__

