#include <set>
#include <stack>
//#include <Eigen/Dense>
#include <armadillo>

#include <gecode/driver.hh>

#include "cbs.hpp"

using namespace Gecode;

using SetView = Set::SetView;
using NaryProp = NaryPropagator<SetView, Set::PC_SET_CGLB>;

class SpanningTreeCtr : public NaryProp {
  using NaryProp::x;
public:
  SpanningTreeCtr(const Home& home, ViewArray<Set::SetView>& x)
    : NaryProp(home, x) {}

  SpanningTreeCtr(Space& home, bool share, SpanningTreeCtr& p)
    : NaryProp(home, share, p) {}

  Actor* copy(Space& home, bool share) override {
    return new (home) SpanningTreeCtr(home,share,*this);
  }

  ExecStatus propagate(Space& home, const ModEventDelta& med) override {
    enforce_node_symmetry(home);

    auto ccs = connected_components_no_cycle<SetVarGlbValues>();
    if (ccs.empty())
      return ES_FAILED;

    prune_LUB_to_avoid_cycle(home, ccs);
    return ES_FIX;
  }

  arma::mat construct_Laplacian_matrix() const {
    int n = x.size();
    arma::mat m(n,n, arma::fill::zeros);
    auto ccs = connected_components_no_cycle<SetVarGlbValues>();

    for (const auto& cc : ccs) {
      auto outDegrees = connected_component_out_degrees(cc);
      int master_node = *cc.begin();
      int master_node_degree = std::accumulate(outDegrees.begin(),
                                               outDegrees.end(), 0);

      for (int i=0; i<n; i++) {
        m(i, master_node) = -outDegrees[i];
        m(master_node, i) = -outDegrees[i];
      }

      m(master_node, master_node) = master_node_degree;

      for (auto node_cc = ++cc.begin(); node_cc != cc.end(); ++node_cc) {
        assert(*node_cc != master_node);
        m(*node_cc, *node_cc) = 1;
      }
    }
    return std::move(m);
  }

  std::vector<int> connected_component_out_degrees(const std::set<int>& cc) const {
    std::vector<int> outDegree;
    outDegree.resize((size_t) x.size());
    std::vector<int> degreeCount_Glb((size_t) x.size(), 0);
    std::vector<int> degreeCount_Lub((size_t) x.size(), 0);
    for (auto node : cc) {
        for (SetVarGlbValues adj(x[node]); adj(); ++adj)
          degreeCount_Glb[adj.val()] += 1;
        for (SetVarLubValues adj(x[node]); adj(); ++adj)
          degreeCount_Lub[adj.val()] += 1;
      }
    for (auto i = 0; i < x.size(); i++)
        outDegree[i] = degreeCount_Lub[i] - degreeCount_Glb[i];
    return std::move(outDegree);
  }

  void slndist(Space& home, SolnDistribution *dist,
               SolnDistribution::Type type) const override {
    arma::mat laplacian = construct_Laplacian_matrix();

    // std::set<n,v> voisins_traitées.
    // FORALL nodes n
      // VAR inverse non computé
      // FORALL voisins v non assignées [SetVarLubValues et notContains]
        // SI n,v non traité
          // SI inverse non computé
            // COMPUTE inverse
          // Assignation densité

    struct Edge {
      int minN; // End node with minimum index
      int maxN; // End node with maximum index
      Edge(int a, int b) {
        minN = std::min(a,b);
        maxN = std::max(a,b);
      }
      bool operator<(const Edge& b) const {
        if (minN == b.minN)
          return maxN < b.maxN;
        return minN < b.minN;
      }
    };

    std::set<Edge> computed;
    for (int i=0; i<x.size(); i++) {
      if (laplacian(i,i) == 1) continue;
      arma::mat inv;
      for_every_incertain_set_values(i, [&](int adj) {
        Edge e{i,adj};
        if (computed.find(e) == computed.end()) {
          computed.insert(e);
          if (inv.empty()) {
            auto idxs = continuous_idx_vector_skip_i(i, x.size());
            inv = arma::inv(laplacian.submat(idxs, idxs));
            std::cout << laplacian << std::endl;
            std::cout << inv << std::endl;
          }
          int jj = i<adj ? adj-1 : adj;
          double dens = inv(jj,jj);
          // C'est affreux, j'envoie seulement une densité sur les deux
          // possibles...
          dist->setMarginalDistribution(id(), x[i].id(), adj, dens);
        }
      });
    }

  }

  void for_every_incertain_set_values(int i,
                                      const std::function<void(int)>& f) const {
    for (SetVarLubValues adj(x[i]); adj(); ++adj) {
      if (!x[i].contains(adj.val())) {
        f(adj.val());
      }
    }

  }

  arma::uvec continuous_idx_vector_skip_i(int i, int n) const {
    arma::uvec idxs((arma::uword)n-1);
    arma::uword num = 0;
    for (int j=0; j<n-1; j++, num++) {
      if (j==i) num++;
      idxs[j] = num;
    }
    return std::move(idxs);
  }


  void
  slndistsize(SolnDistributionSize *s, unsigned int& domSum,
                   unsigned int& domSumB) const override {
    domSum = 0;
    domSumB = 0;
    for (int i=0; i<x.size(); i++) {
      if (!x[i].assigned()) {
        domSum += x[i].lubSize() - x[i].glbSize();
        if (s->varInBrancher(x[i].id()))
          domSumB += x[i].lubSize() - x[i].glbSize();
      }
    }
  }

  static ExecStatus post(Home home, ViewArray<Set::SetView>& x) {
    (void) new (home) SpanningTreeCtr(home, x);
    return ES_OK;
  }
private:
  template<class SetVarValues>
  std::vector<std::set<int>>
  connected_components_no_cycle() const {
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
          if (adj.val() != parent_n)
            s.push_back({n,adj.val()});
      }
    }
    return std::move(connected_components);
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
  std::set<int>
  set_of_all_nodes(int n) const {
    std::set<int> not_visited;
    for (int i = 0; i < n; i++)
      not_visited.insert(not_visited.end(), i);
    return std::move(not_visited);
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

//    branch(*this, nodes, SET_VAR_SIZE_MIN(), SET_VAL_MIN_INC());
    cbsbranch(*this, nodes,  CBSBranchingHeuristic::MAX_SD);
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
