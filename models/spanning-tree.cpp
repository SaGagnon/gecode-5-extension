#include <set>
#include <stack>
#include <armadillo>

#include <gecode/driver.hh>
#include <unordered_map>

#include "cbs.hpp"

/**
 * CONFIG VARIABLES
 */
bool ALGO_APPROX = false;
bool LAPLACIAN_TRIM = false;

// Pour voir comment l'algo de calcul exact réagit à une imprécision.
//#define PERTURBATION

#ifdef PERTURBATION
/**
 * RANDOM
 */
std::mt19937 gen;
#endif

// Pour voir les densités moyenne calculées
//#define MOY

#ifdef MOY
double _DENS_TOT = 0;
double _CNT = 0;
#endif

using namespace Gecode;

using node_t      = unsigned int;
using edge_t      = std::pair<node_t, node_t>;

using n_nodes_t   = unsigned int;
using n_edges_t   = unsigned int;

namespace io {
  std::vector<std::string> lines(const std::string& file) {
    std::vector<std::string> lines;
    std::ifstream f(file);
    std::string line;
    while (getline(f, line)) lines.push_back(line);
    return std::move(lines);
  }

  std::tuple<n_nodes_t, n_edges_t, std::vector<edge_t>>
  graph(const std::vector<std::string>& lines) {
    const bool TP = false;
    assert(lines.size() >= 3);
    auto n_nodes = stoi(lines[0]);
    auto n_edges = stoi(lines[1]);

    // When reading graphs from TP3 INF4705 H17
    int excess = 0;
    if (lines.size() == 2 + n_nodes + n_edges)
      excess = n_nodes;

    assert(lines.size() == 2 + excess + n_edges);

    std::vector<edge_t> edges((size_t)n_edges);
    transform(lines.begin() + 2 + excess, lines.end(), edges.begin(),
              [](const std::string& line) {
                int n1, n2;
                std::stringstream(line) >> n1 >> n2;
                return std::make_pair(n1 - 1, n2 - 1);
              });

    return std::make_tuple(n_nodes, n_edges, std::move(edges));
  };
}

namespace utils {
  /**
   * Construct set of T from 0 to n
   */
  template<class T>
  std::set<T> set_zero_to(unsigned int n) {
    std::set<T> not_visited;
    for (T i=0; i<n; i++)
      not_visited.insert(not_visited.end(), i);
    return std::move(not_visited);
  }
}

using BoolView = Int::BoolView;
using NaryProp = NaryPropagator<BoolView, Int::PC_BOOL_VAL>;

class SpanningTreeCtr : public NaryProp {
protected:
  using node_graph_t = std::pair<node_t, BoolView>;
  // Adjacency list
  struct {
    node_graph_t** x; // List of adjacent nodes for each node
    unsigned int* size; // Size of adj. list for each node
    n_nodes_t n_nodes;
  } graph{nullptr, nullptr, 0};
public:
  SpanningTreeCtr(Space& home, ViewArray<BoolView>& e,
                  const std::vector<edge_t>& edges,
                  n_nodes_t n_nodes)
    : NaryProp(home, e) {
    graph.n_nodes = n_nodes;
    graph.x = home.alloc<node_graph_t*>(n_nodes);

    auto iter_edge =
      [&](std::function<void(BoolView&, node_t, node_t)> f) {
      for (unsigned int i=0; i<edges.size(); i++) {
        node_t n1, n2;
        std::tie(n1, n2) = edges[i];
        assert(n1 >= 0 && n1 < n_nodes);
        assert(n2 >= 0 && n2 < n_nodes);
        f(e[i], n1, n2);
      }
    };

    graph.size = home.alloc<unsigned int>(n_nodes);
    for (int i=0; i<n_nodes; i++) {
      assert(graph.size[i] == 0);
    }

    iter_edge([&](BoolView&, node_t n1, node_t n2) {
      graph.size[n1] += 1;
      graph.size[n2] += 1;
    });

    for (int i=0; i<n_nodes; i++) {
      graph.x[i] = home.alloc<node_graph_t>(graph.size[i]);
    }

    Region r(home);
    auto ins_pos = r.alloc<unsigned int>(n_nodes);
    for (int i=0; i<n_nodes; i++) {
      assert(ins_pos[i] == 0);
    }

    iter_edge([&](BoolView& view, node_t n1, node_t n2) {
      graph.x[n1][ins_pos[n1]++] = std::make_pair(n2, view);
      graph.x[n2][ins_pos[n2]++] = std::make_pair(n1, view);
    });

    for (int i=0; i<n_nodes; i++) {
      assert(ins_pos[i] == graph.size[i]);
    }
  }

  SpanningTreeCtr(Space& home, bool share, SpanningTreeCtr& p)
    : NaryProp(home, share, p) {
    graph.n_nodes = p.graph.n_nodes;
    graph.size    = home.alloc<unsigned int> (graph.n_nodes);
    graph.x       = home.alloc<node_graph_t*>(graph.n_nodes);

    for (int i=0; i<graph.n_nodes; i++) {
      graph.size[i] = p.graph.size[i];
    }

    for (int i=0; i<graph.n_nodes; i++) {
      graph.x[i] = home.alloc<node_graph_t>(graph.size[i]);
      for (int j=0; j<graph.size[i]; j++) {
        graph.x[i][j].first = p.graph.x[i][j].first;
        graph.x[i][j].second.update(home, share, p.graph.x[i][j].second);
      }
    }
  }

  virtual PropCost cost(const Space& home, const ModEventDelta& med) const {
    return PropCost::cubic(PropCost::HI,graph.n_nodes);
  }

  static ExecStatus post(Home home, ViewArray<BoolView>& e,
                         const std::vector<edge_t>& edges,
                         n_nodes_t n_nodes) {
    auto p = new (home) SpanningTreeCtr(home, e, edges, n_nodes);
    // The propagator is not scheduled for propagation at the beginning, so
    // we have to do it by hands (because of type PC_BOOL_VAL). We can't change
    // for PC_INT_DOM, as it gets automatically mapped to PC_BOOL_VAL when
    // the propagator subscribes to its variables during creation.
    // (see section 23.1 of doc). Value 0 of second argument is dummy.
    return  p->propagate(home, 0);
  }

  Actor* copy(Space& home, bool share) override {
    return new (home) SpanningTreeCtr(home,share,*this);
  }

  ExecStatus propagate(Space& home, const ModEventDelta& med) override {

    // First, all nodes must be accessible via unnassigned edges and
    // edges equal to one. If this is not the case, we can't have
    // a connected graph
    {
      std::set<node_t> seen_nodes;
      DFS_all(0, [&](node_t node) {
        seen_nodes.insert(node);
        return ES_OK;
      });
      if (seen_nodes.size() != graph.n_nodes) {
        return ES_FAILED;
      }
    }

    // unvisited nodes
    auto u_nodes = utils::set_zero_to<node_t>(graph.n_nodes);
    while (!u_nodes.empty()) {
      std::set<node_t> CC;

      // Construction of CC (connected component) without cycles
      auto ret = DFS_one(*u_nodes.begin(), [&](node_t node) {
        // If this fails, we have a cycle.
        if (!CC.insert(node).second)
          return ES_FAILED;
        u_nodes.erase(node);
        return ES_OK;
      });
      if (ret != ES_OK) return ret;


      if (CC.size() >= 3) {
        // We make sure there's no unassigned edges pointing to the CC itself
        ret = DFS_one(*CC.begin(), [&](node_t node) {
          for_each_adj(node, [&](node_t adj_n, BoolView& adj_v) {
            if (adj_v.none() && CC.find(adj_n) != CC.end()) {
              adj_v.eq(home, 0);
              // See last comment
              return ES_NOFIX;
            }
          });
          return ES_OK;
        });
        if (ret != ES_OK) return ret;
      }


      // Number of adjacent unassigned edges for the cc
      int n_unassigned = 0;
      BoolView last_view_unassigned;

      ret = DFS_one(*CC.begin(), [&](node_t node) {
        for_each_adj(node, [&](node_t adj_n, BoolView& adj_v) {
          if (adj_v.none()) {
            n_unassigned++;
            last_view_unassigned = adj_v;
          }
        });
        return ES_OK;
      });
      if (ret != ES_OK) return ret;

      // If there's no adjacent edges, we fail
      if (n_unassigned == 0 && CC.size() != graph.n_nodes)
        return ES_FAILED;
      // If there's only one, we are obliged to take it
      if (n_unassigned == 1) {
        last_view_unassigned.eq(home, 1);
        // If we fix a variable, we have to redo the propagation process;
        // CCs might change and trigger other changes. This is dirty and
        // probably very inneficient. However, I'm only interested in getting
        // propagation to work correctly right now. I'm sure it can be
        // optimized.
        return ES_NOFIX;
      }

    }

    if (x.assigned()) {
      return home.ES_SUBSUMED(*this);
    }
    return ES_FIX;
  }

  void solndistrib(Space& home, SolnDistrib *dist) const override {

    #ifdef PERTURBATION
    std::uniform_real_distribution<> dis(-0.3, 0);
//    std::uniform_real_distribution<> dis(-0.05, 0.05);
    #endif

    // each node belong to a connected component identified by a master node
    std::unordered_map<node_t, unsigned int> ncc;
    {
      // unvisited nodes
      auto u_nodes = utils::set_zero_to<node_t>(graph.n_nodes);
      while (!u_nodes.empty()) {
        auto m_node = *u_nodes.begin();
        DFS_one(m_node, [&](node_t node) {
          u_nodes.erase(node);
          ncc[node] = m_node;
          return ES_OK;
        });
      }
    }


    arma::mat laplacian;
    {
      std::set<int> ___seen_nodes;
      laplacian = arma::mat(graph.n_nodes, graph.n_nodes, arma::fill::zeros);
      DFS_all(0, [&](node_t node) {
        assert(___seen_nodes.find(node) == ___seen_nodes.end());
        ___seen_nodes.insert(node);
        if (ncc[node] != node)
          laplacian(node, node) = 1;
        for_each_adj(node, [&](node_t adj_n, BoolView& adj_v) {
          // We are going to see each edge two times and we only want
          // to add it one time
          if (!adj_v.none()) return ES_OK;
          if (adj_v.none() && node < adj_n) {
            auto cc1 = ncc[node];
            auto cc2 = ncc[adj_n];
            assert(cc1 != cc2);
            laplacian(cc1, cc2) -= 1;
            laplacian(cc2, cc1) -= 1;
            laplacian(cc1, cc1) += 1;
            laplacian(cc2, cc2) += 1;
          }
        });
        return ES_OK;
      });
    }

    if (!ALGO_APPROX) {
      Region r(home);
      // Number of edges whose density is not calculated per node
      auto *cards = r.alloc<unsigned int>(graph.n_nodes);
      for (int jj = 0; jj < graph.n_nodes; jj++) {
        unsigned int degree = laplacian(jj, jj);
        if (degree == 1) { // Because of contraction
          cards[jj] = 0;
        } else {
          cards[jj] = degree;
        }
      }

      // CC refering to master nodes to their position in the reduced laplacian
      std::unordered_map<unsigned, unsigned> cc_to_minilap;
      {
        unsigned pos = 0;
        for(int i=0; i<graph.n_nodes; i++)
          if (cards[i] != 0)
            cc_to_minilap[i] = pos++;

        std::cout << "largeur = " << pos << std::endl;
      }


      // Visited edge ids (for which we already calculated densities)
      std::set<unsigned int> visitied_edges;

      while (true) {
        // Selected jj for inverting matrix and assigning densities
        unsigned int jj;
        {
          auto nptr = std::max_element(cards, cards + graph.n_nodes);
          if (*nptr == 0) break;
          jj = (node_t) (nptr - cards);
        }

        arma::mat inv;
        {
          arma::uvec idxs;
          if (LAPLACIAN_TRIM) {
            std::vector<unsigned> _idxs;
            unsigned n = 0;
            for (int i=0; i<graph.n_nodes; i++)
              if (i != jj && laplacian(i,i) != 1)
                _idxs.push_back(i);

            // TODO: REFAIRE ÇA...
            idxs = arma::uvec(_idxs.size());
            for (int i=0; i<_idxs.size(); i++) {
              idxs[i] = _idxs[i];
            }
          } else {
            idxs = arma::uvec(graph.n_nodes - 1);
            int n = 0;
            for (int i = 0; i < graph.n_nodes - 1; i++) {
              if (i == jj) n++;
              idxs[i] = n++;
            }
          }
          // TODO: Regarder si on a un slow down a cause du view
          inv = arma::inv_sympd(laplacian.submat(idxs, idxs));
        }

        // We have to explore each edge in each node of the cc jj
        for (int n = 0; n < graph.n_nodes; n++) {
          if (ncc[n] == jj) {
            for_each_adj(n, [&](node_t adj_n, BoolView& adj_v) {
              if (adj_v.none() &&
                  visitied_edges.find(adj_v.id()) == visitied_edges.end()) {
                // Diagonal index for density
                auto ii = ncc[adj_n];
                // Because size(laplacian) > size(inv)
                double dens;

                assert(cc_to_minilap.find(ii) != cc_to_minilap.end());

                // ii and jj will change depending on if we reduce the laplacian
                // or not when inverting.
                unsigned lii, ljj;
                if (LAPLACIAN_TRIM) {
                  lii = cc_to_minilap[ii];
                  ljj = cc_to_minilap[jj];
                } else {
                  lii = ii;
                  ljj = jj;
                }

                 if (lii > ljj) {
                  dens = inv(lii-1,lii-1);
                } else {
                  dens = inv(lii, lii);
                }

                #ifdef PERTURBATION
                dens += dis(gen);
                #endif

                #ifdef MOY
                _DENS_TOT += dens;
                _CNT++;
                #endif

                dist->marginaldistrib(id(), adj_v.id(), 1, dens);
                dist->marginaldistrib(id(), adj_v.id(), 0, 1 - dens);
                visitied_edges.insert(adj_v.id());

                assert(cards[ii] !=0);
                assert(cards[jj] !=0);
                cards[jj] -= 1;
                cards[ii] -= 1;
              }
            });
          }
        }
      }

    } else {

      const int n = graph.n_nodes;

      // Compute f[L,i], forall nodes
      Region r(home);
      auto fL = r.alloc<int>(n);
      {
        for (int i = 0; i < n; i++) {
          assert(fL[i] == 0);
        }

        for (int i = 0; i < n; i++) {
          for (int k = 0; k < n - 1; k++)
            if (k != i)
              for (int l = k + 1; l < n; l++)
                if (l != i)
                  fL[i] += 2 * laplacian(k,i) * laplacian(k,l) * laplacian(l,i);
          for (int k = 0; k < n; k++)
            if (k != i)
              fL[i] += std::pow(laplacian(k,i), 2) * laplacian(k, k);
        }
      }

      // We visit each edge and assign its density
      DFS_all(0, [&](node_t node) {
        auto jj = ncc[node];
        for_each_adj(node, [&](node_t adj_n, BoolView& adj_v) {
          if (adj_v.none()) {
            auto ii = ncc[adj_n];
            // We are going to see each edge two times and we only want
            // to add it one time
            if (ii > jj) {
              /**
               * Body itération chaque edge (i,j)
               */
              // Called f[A,i], maybe change name in article.
              auto A = fL[ii];
              // On enlève tous les occurences lorsque l == j, forall k
              for (int k = 0; k < n; k++)
                if (k != ii)
                  A -= laplacian(k, ii) * laplacian(k, jj) * laplacian(jj, ii);
              // On enlève tous les occurences lorsque k == j, forall l
              for (int l = 0; l < n; l++)
                if (l != ii)
                  A -= laplacian(jj, ii) * laplacian(jj, l) * laplacian(l, ii);

              // On a un élément deux fois
              A += laplacian(ii, jj) * laplacian(jj, jj) * laplacian(jj, ii);

              int B = 0;
              for (int k = 0; k < n; k++)
                if (k != ii && k != jj)
                  B += std::pow(laplacian(k, ii), 2);

              double dens;
              if (A == 0)
                dens = 1 / laplacian(ii,ii);
              else
                dens = (double) A / (laplacian(ii, ii) * A - std::pow(B, 2));

              #ifdef PERTURBATION
              dens += dis(gen); // TEST -- TMP
              #endif

              #ifdef MOY
              _DENS_TOT += dens;
              _CNT++;
              #endif

              dist->marginaldistrib(id(), adj_v.id(), 1, dens);
              dist->marginaldistrib(id(), adj_v.id(), 0, 1 - dens);
            }
          }
        });
        return ES_OK;
      });

    }
  }

  void solndistribsize(SolnDistribSize *s, unsigned int& domsum,
              unsigned int& domsum_b) const override {
    domsum = 0;
    domsum_b = 0;
    for (auto var : x) {
      if (!var.assigned()) {
        domsum += var.size();
        if (s->inbrancher(var.id()))
          domsum_b += var.size();
      }
    }
  }
private:
  void for_each_adj(node_t node,
                    const std::function<void(node_t,BoolView&)>& f) const {
    for (int i=0; i<graph.size[node]; i++) {
      node_t adj_node; BoolView adj_view;
      std::tie(adj_node, adj_view) = graph.x[node][i];
      f(adj_node, adj_view);
    }
  }

  /**
   * Depth first search
   *
   * DFS in graph that considers only assigned edges in graph
   *
   * @param start starting node
   * @param visit function to apply to each node
   * @return
   */
  ExecStatus DFS_one(node_t start,
                     const std::function<ExecStatus(node_t)>& visit) const {
    using parent_node_t = node_t;
    std::stack<std::pair<parent_node_t, node_t>> stack;
    stack.emplace(start, start);
    while (!stack.empty()) {
      parent_node_t parent_node; node_t node;
      std::tie(parent_node, node) = stack.top();
      stack.pop();
      auto ret = visit(node);
      if (ret != ES_OK) return ret;
      for_each_adj(node, [&](node_t adj_node, BoolView& adj_view) {
        if (adj_view.one() && adj_node != parent_node)
          stack.emplace(node, adj_node);
      });
    }
    return ES_OK;
  }

  /**
   * Depth first search
   *
   * DFS considering all assigned and unassigned edges in graph.
   *
   * @param start see above
   * @param visit see above
   * @return
   */
  ExecStatus DFS_all(node_t start,
                     const std::function<ExecStatus(node_t)>& visit) const {
    std::set<node_t> visited;
    std::stack<node_t> stack;
    stack.push(start);
    while (!stack.empty()) {
      auto node = stack.top();
      stack.pop();
      // If we did not already visit this node
      if (visited.insert(node).second) {
        auto ret = visit(node);
        if (ret != ES_OK) return ret;
        for_each_adj(node, [&](node_t adj_node, BoolView& adj_view) {
          if ((adj_view.one() || adj_view.none())
              && visited.find(adj_node) == visited.end()) {
            stack.push(adj_node);
          }
        });
      }
    }
    return ES_OK;
  }

};

void spanning_tree_ctr(Space& home, const BoolVarArgs& e,
                       const std::vector<edge_t>& edges,
                       n_nodes_t n_nodes) {
  ViewArray<BoolView> _e(home, e);

  rel(home, sum(e) == n_nodes - 1);

  SpanningTreeCtr::post(home, _e, edges, n_nodes);
}

class ConstrainedSpanningTree : public Script {
private:
protected:
  BoolVarArray e;
  const std::string instance;
public:
  explicit ConstrainedSpanningTree(const InstanceOptions& opt)
    : Script(opt), instance(opt.instance()) {

    n_nodes_t           n_nodes;
    n_edges_t           n_edges;
    std::vector<edge_t> edges;

    std::tie(n_nodes, n_edges, edges) = io::graph(io::lines(opt.instance()));

    e =  BoolVarArray(*this, (int)edges.size(), 0, 1);

    spanning_tree_ctr(*this, e, edges, n_nodes);

    std::vector<BoolVarArgs> adj_edges(n_nodes);
    assert(adj_edges.size() == n_nodes);
    for (int i=0; i<n_edges; i++) {
      node_t n1, n2;
      std::tie(n1, n2) = edges[i];

      adj_edges[n1] << e[i];
      adj_edges[n2] << e[i];
    }

    for (int i=0; i<n_nodes; i++) {
      assert(adj_edges[i].size() != 0);
      rel(*this, sum(adj_edges[i]) <= 2);
    }

    cbsbranch(*this, e,  CBSBranchingHeuristic::MAX_SD);
    branch(*this, e, BOOL_VAR_NONE(), BOOL_VAL_MIN());
  }

  ConstrainedSpanningTree(bool share, ConstrainedSpanningTree& s)
    : Script(share,s), instance(s.instance) {
    e.update(*this, share, s.e);
  }

  Space* copy(bool share) override {
    return new ConstrainedSpanningTree(share,*this);
  }
};

int
main(int argc, char* argv[]) {
  #ifdef PERTURBATION
  std::random_device rd;
  gen = std::mt19937(rd());
  #endif

  InstanceOptions opt("Spanning Tree");
  opt.ipl(IPL_DOM);
  opt.solutions(1);
//  opt.mode(SM_GIST);
  opt.parse(argc, argv);

  for (int i=0; i<argc; i++) {
    std::string curr_param = argv[i];
    if (curr_param == "-approx")
      ALGO_APPROX = true;
    else if (curr_param == "-trim")
      LAPLACIAN_TRIM = true;
  }
  Script::run<ConstrainedSpanningTree,DFS,InstanceOptions>(opt);

  #ifdef MOY
  std::cout << "DENS MOY = " << _DENS_TOT / _CNT << std::endl;
  #endif


}
