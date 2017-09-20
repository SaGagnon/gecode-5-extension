#include <set>

#include <gecode/driver.hh>
#include <stack>

//#include "cbs.hpp"

using namespace Gecode;

using SetView = Set::SetView;
using NaryProp = NaryPropagator<SetView, Set::PC_SET_CGLB>;

class SpanningTreeCtr : public NaryProp {
  using NaryProp::x;
private:
  template<class SetVarValues>
  std::vector<std::set<int>>
  connected_components_no_cycle() {
    using CC = std::vector<std::set<int>>;
    using Node = int;
    using SetNode = std::set<int>;
    using ParentNode_Node = std::pair<int,int>;

    CC connected_components;
    auto all_nodes = set_of_all_nodes(x.size());
    // DFS FORALL connected components
    while (!all_nodes.empty()) {
      // We begin DFS with first_node
      auto first_node = all_nodes.begin();
      std::deque<ParentNode_Node > s{{-1,*first_node}};
      connected_components.emplace_back(SetNode{});
      // DFS
      while (!s.empty()) {
        // pop from stack
        auto pair = s.back(); s.pop_back();
        Node parent_n = pair.first;
        Node n = pair.second;
        // we have a cycle if n is not in all_nodes
        if (all_nodes.find(n) != all_nodes.end()) {
          all_nodes.erase(n);
          connected_components.back().insert(n);
        } else {
          return CC{};
        }
        // push childs except parent
        for (SetVarValues adj(x[n]); adj(); ++adj)
          if (adj.val() != parent_n && adj.val() != n)
            s.push_back({n,adj.val()});
      }
    }
    return std::move(connected_components);
  }

public:
  SpanningTreeCtr(const Home& home, ViewArray<Set::SetView>& x)
    : NaryProp(home, x) {}

  SpanningTreeCtr(Space& home, bool share, SpanningTreeCtr& p)
    : NaryProp(home, share, p) {}

  Actor* copy(Space& home, bool share) override {
    return new (home) SpanningTreeCtr(home,share,*this);
  }

  void prune_LUB_to_avoid_cycle(Space& home, std::vector<std::set<int>> ccs) {
    for (auto cc : ccs)
      for (auto node : cc)
        for (auto node_to_delete : cc)
          if (!x[node].contains(node_to_delete))
            x[node].exclude(home, node_to_delete);
  }

  void enforce_node_symmetry(Space& home) {
    for (int i=0; i<x.size(); i++) {
      for (SetVarGlbValues adj(x[i]); adj(); ++adj) {
        x[adj.val()].include(home, i);
      }
    }
  }

  ExecStatus propagate(Space& home, const ModEventDelta& med) override {
    enforce_node_symmetry(home);

    auto ccs = connected_components_no_cycle<SetVarGlbValues>();
    if (ccs.empty())
      return ES_FAILED;

    prune_LUB_to_avoid_cycle(home, ccs);
    return ES_FIX;
  }

  std::set<int>
  set_of_all_nodes(int n) const {
    std::set<int> not_visited;
    for (int i = 0; i < n; i++)
      not_visited.insert(not_visited.end(), i);
    return std::move(not_visited);
  }

  void slndist(Space& home, SolnDistribution *dist,
               SolnDistribution::Type type) const override {
    Propagator::slndist(home, dist, type);
  }
  void
  slndistsize(SolnDistributionSize *s, unsigned int& domAggr,
                   unsigned int& domAggrB) const override {
    Propagator::slndistsize(s, domAggr, domAggrB);
  }

  static ExecStatus post(Home home, ViewArray<Set::SetView>& x) {
    (void) new (home) SpanningTreeCtr(home, x);
    return ES_OK;
  }
};

void spanningTreeCtr(Home home, const SetVarArgs& x) {
  int n = x.size();

  IntVarArgs cards(home, n, 1, 2);
  for (int i=0; i<n; i++)
    rel(home, cards[i] == cardinality(x[i]));
  rel(home, sum(cards)/2 == n - 1 );
  channel(home, x, x);

  ViewArray<SetView> y(home,x);
  GECODE_ES_FAIL(SpanningTreeCtr::post(home,y));
}

class SpanningTree : public Script {
private:
  // UTILITIES
  using Edge = std::pair<int, int>;
  std::vector<std::string> read_file_lines(const InstanceOptions& opt) const {
    std::vector<std::string> lines;
    std::ifstream file(opt.instance());
    std::string line;
    while (getline(file, line)) lines.push_back(line);
    return lines;
  }
  std::vector<std::vector<int>>
  adj_matrix(int n_nodes, const std::vector<Edge>& edges) const {
    std::vector<std::vector<int>> adj((size_t)n_nodes);
    for (auto edge : edges) {
      int n1 = edge.first, n2 = edge.second;
      adj[n1].push_back(n2);
      adj[n2].push_back(n1);
    }
    return adj;
  }
  std::vector<Edge>
  get_edges_from_file(std::vector<std::string>& lines,
                      int n_nodes, int n_edges) const {
    std::vector<std::pair<int,int>> edges((size_t)n_edges);
    transform(lines.begin()+2+n_nodes, lines.end(), edges.begin(),
              [](std::string line) {
                int n1, n2; std::stringstream(line) >> n1 >> n2;
                return std::make_pair(n1 - 1, n2 - 1);
              });
    return edges;
  }
protected:
  SetVarArray nodes;
public:
  explicit SpanningTree(const InstanceOptions& opt)
    : Script(opt) {
    auto lines = read_file_lines(opt);

    int n_nodes = stoi(lines[0]);
    int n_edges = stoi(lines[1]);

    auto edges = get_edges_from_file(lines, n_nodes, n_edges);
    auto adj = adj_matrix(n_nodes, edges);

    // nodes[i] contient tous ses noeuds adjacents.
    nodes = SetVarArray(*this, n_nodes, IntSet::empty, IntSet(0, n_nodes-1), 1, 2);
    for(int i=0; i<adj.size(); i++)
      dom(*this, nodes[i], SRT_SUB, IntSet(adj[i].data(), (int)adj[i].size()));

    spanningTreeCtr(*this, nodes);

    branch(*this, nodes, SET_VAR_SIZE_MIN(), SET_VAL_MIN_INC());
  }

  SpanningTree(bool share, SpanningTree& s)
    : Script(share,s) {
    nodes.update(*this, share, s.nodes);
  }

  Space* copy(bool share) override {
    return new SpanningTree(share,*this);
  }

  void print(std::ostream& os) const override {
    for (int i=0; i<nodes.size(); i++) {
     os << nodes[i] << std::endl;
    }
  }
};

int
main(int argc, char* argv[]) {
  InstanceOptions opt("Spanning Tree");
  opt.ipl(IPL_DOM);
  opt.solutions(1);
//  opt.mode(SM_GIST);

  opt.parse(argc, argv);

  Script::run<SpanningTree,DFS,InstanceOptions>(opt);
}
