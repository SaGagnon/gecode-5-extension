#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/minimodel.hh>

#include <fstream>
#include <sstream>
#include <string>

#include <cbs.hpp>

using namespace Gecode;

class MinDef : public IntMinimizeScript {
protected:
  std::string file_path;

  int num_vertices;
  int num_edges;

//  static int *cardinalities;

  IntVarArray edge_col;
  IntVar n_holes;
public:
  MinDef(const InstanceOptions& opt)
    : IntMinimizeScript(opt), file_path(opt.instance()) {

    std::vector<std::vector<int>> _node_edges((unsigned long) num_vertices);
    {
      std::ifstream file(file_path);
      file >> num_vertices;
      file >> num_edges;

      // Go to first line of data
      {
        std::string tmp;
        std::getline(file, tmp);
      }

      _node_edges.resize((unsigned long) num_vertices);
      for (int i = 0; i < num_vertices; i++) {
        std::string line;
        std::getline(file, line);
        std::istringstream iss(line);

        int edge;
        while (iss >> edge)
          _node_edges[i].push_back(edge);
      }
    }

    int n_cols;
    {
      int max_card = 0;
      for (int i=0; i<num_vertices; i++) {
        auto card = (int)_node_edges[i].size();
        if (card > max_card)
          max_card = card;
      }

      n_cols = max_card*2;
      edge_col = IntVarArray(*this, num_edges, 0, max_card);
    }


    IntVarArgs node_holes(*this, num_vertices, 0, n_cols-2);
    for (int i=0; i<num_vertices; i++) {
      auto card = (int)_node_edges[i].size();
      IntVarArgs v_edges(card);
      for (int j=0; j<card; j++)
        v_edges[j] = edge_col[_node_edges[i][j]-1];
      distinct(*this, v_edges, opt.ipl());
      rel(*this, node_holes[i] == max(v_edges) - min(v_edges) + 1 - card);
    }

    n_holes = IntVar(*this, 0, num_edges);
    rel(*this, n_holes == sum(node_holes));

    cbsbranch(*this, edge_col, CBSBranchingHeuristic::MAX_SD);
//    branch(*this, edge_col, INT_VAR_SIZE_MIN(), INT_VAL_SPLIT_MIN());
//    branch(*this, edge_col, INT_VAR_AFC_MIN(opt.decay()), INT_VAL_SPLIT_MIN());
//    branch(*this, edge_col, INT_VAR_ACTIVITY_SIZE_MAX(opt.decay()), INT_VAL_MIN());
  }

  MinDef(bool share, MinDef& s)
    : IntMinimizeScript(share,s) {
    edge_col.update(*this,share,s.edge_col);
    n_holes.update(*this,share,s.n_holes);
  }

  virtual Space*
  copy(bool share) {
    return new MinDef(share,*this);
  }

  virtual IntVar cost(void) const {
    return n_holes;
  }

  virtual void
  print(std::ostream& os) const {
    os << "Score: " << n_holes << std::endl;
  }
};

int main(int argc, char *argv[]) {
  InstanceOptions opt("Minimum Deficiency");
  opt.ipl(IPL_DOM);
  opt.solutions(0);

  opt.parse(argc,argv);

  Script::run<MinDef,BAB,InstanceOptions>(opt);
  return 0;
}
