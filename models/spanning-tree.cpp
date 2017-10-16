#include <set>
#include <stack>
#include <armadillo>

#include <gecode/driver.hh>

#include "cbs.hpp"

using namespace Gecode;

using node_t      = unsigned int;
using edge_t      = std::pair<node_t, node_t>;

using n_nodes_t   = unsigned int;
using n_edges_t   = unsigned int;

//using adj_list_t  = std::vector<std::vector<node_t>>;

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
    assert(lines.size() >= 3);
    auto n_nodes = stoi(lines[0]);
    auto n_edges = stoi(lines[1]);

    assert(lines.size() == 2 + n_nodes + n_edges);

    std::vector<edge_t> edges((size_t)n_edges);
    transform(lines.begin() + 2 + n_nodes, lines.end(), edges.begin(),
              [](const std::string& line) {
                int n1, n2;
                std::stringstream(line) >> n1 >> n2;
                return std::make_pair(n1 - 1, n2 - 1);
              });

    return std::make_tuple(n_nodes, n_edges, std::move(edges));
  };

//  adj_list_t
//  adj_list(n_nodes_t n_nodes, const std::vector<edge_t>& edges) {
//    assert(n_nodes > 1);
//    assert(edges.size() > 1);
//    adj_list_t adj(n_nodes);
//
//    for (auto edge : edges) {
//      node_t n1, n2;
//      std::tie(n1, n2) = edge;
//
//      adj[n1].push_back(n2);
//      adj[n2].push_back(n1);
//    }
//
//    for (auto list : adj) {
//      assert(list.size() > 0);
//    }
//
//    return std::move(adj);
//  }
}

namespace utils {
  template<class T>
  std::set<T> set_zero_to(unsigned int n) {
    std::set<T> not_visited;
    for (T i=0; i<n; i++)
      not_visited.insert(not_visited.end(), i);
    return std::move(not_visited);
  }
}

using BoolView = Int::BoolView;
using NaryProp = NaryPropagator<BoolView, Int::ME_INT_DOM>;

class SpanningTreeCtr : public NaryProp {
protected:
  using node_graph_t = std::pair<node_t, BoolView>;
  // Adjacency list
  struct {
    node_graph_t** x; // List of adjacent nodes for each node
    unsigned int* size; // Size of adj. list for each node
    n_nodes_t n_nodes;
  } graph;
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

  static ExecStatus post(Home home, ViewArray<BoolView>& e,
                         const std::vector<edge_t>& edges,
                         n_nodes_t n_nodes) {
    (void) new (home) SpanningTreeCtr(home, e, edges, n_nodes);
    return ES_OK;
  }

  Actor* copy(Space& home, bool share) override {
    return new (home) SpanningTreeCtr(home,share,*this);
  }

  ExecStatus propagate(Space& home, const ModEventDelta& med) override {

    using visit_func = std::function<ExecStatus(node_t)>;
    using adj_func = std::function<ExecStatus(node_t,BoolView&)>;
    auto DFS = [&](node_t start, visit_func visit, adj_func adj) {
      using parent_node_t = node_t;
      std::deque<std::pair<parent_node_t, node_t>> stack;
      stack.emplace_back(start, start);
      while (!stack.empty()) {
        parent_node_t parent_node; node_t node;
        std::tie(parent_node, node) = stack.back();
        stack.pop_back();
        auto ret = visit(node);
        if (ret != ES_OK) return ret;
        for (int i=0; i<graph.size[node]; i++) {
          node_t adj_node; BoolView adj_view;
          std::tie(adj_node, adj_view) = graph.x[node][i];
          if (adj_view.one() && adj_node != parent_node) {
            stack.emplace_back(node, adj_node);
          }
          ret = adj(adj_node, adj_view);
          if (ret != ES_OK) return ret;
        }
      }
    };


    // unvisited nodes
    auto u_nodes = utils::set_zero_to<node_t>(graph.n_nodes);
    while (!u_nodes.empty()) {
      std::set<node_t> CC;

      // Construction of CC (connected component) without cycles
      auto ret = DFS( *u_nodes.begin(),
                      [&](node_t node) {
                        if (!CC.insert(node).second)
                          return ES_FAILED;
                        u_nodes.erase(node);
                        return ES_OK;
                      },
                      [](node_t,BoolView&) { return ES_OK; }
      );
      if (ret != ES_OK) return ret;

      // If there's only one node
      if (CC.size() == 1) {

        // Number of adjacent unassigned edges
        int n_unassigned = 0;
        BoolView last_view_unassigned;

        ret = DFS( *CC.begin(),
                   [](node_t) { return ES_OK; },
                   [&](node_t, BoolView& adj_v) {
                     if (adj_v.none()) {
                       n_unassigned++;
                       last_view_unassigned = adj_v;
                     }
                     return ES_OK;
                   }
        );
        if (ret != ES_OK) return ret;

        // If there's no adjacent edges, we fail
        if (n_unassigned == 0)
          return ES_FAILED;
        // If there's only one, we are obliged to take it
        if (n_unassigned == 1)
          last_view_unassigned.eq(home, 1);

      } else if (CC.size() > 3) {

        // We make sure there's no unassigned edges pointing to the CC itself
        ret = DFS(*CC.begin(),
                  [](node_t) { return ES_OK; },
                  [&](node_t adj_n, BoolView& adj_v) {
                    if (adj_v.none() && CC.find(adj_n) != CC.end())
                      adj_v.eq(home, 0);
                    return ES_OK;
                  }
        );
        if (ret != ES_OK) return ret;
      }
    }

    return ES_FIX;
  }

  void slndist(Space& home, SlnDist *dist) const override {
  }

  void slndistsize(SlnDistSize *s, unsigned int& domSum,
              unsigned int& domSumB) const override {
    domSum = 0;
    domSumB = 0;
    for (int i=0; i<x.size(); i++) {
      if (!x[i].assigned()) {
        domSum += x[i].size();
        if (s->inbrancher(x[i].id()))
          domSumB += x[i].size();
      }
    }
  }
};

void spanning_tree_ctr(Space& home, const BoolVarArgs& e,
                       const std::vector<edge_t>& edges,
                       n_nodes_t n_nodes) {
  ViewArray<BoolView> _e(home, e);
  SpanningTreeCtr::post(home, _e, edges, n_nodes);
}

class SpanningTree : public Script {
private:
protected:
  BoolVarArray e;
public:
  explicit SpanningTree(const InstanceOptions& opt)
    : Script(opt) {

    n_nodes_t           n_nodes;
    n_edges_t           n_edges;
    std::vector<edge_t> edges;

    std::tie(n_nodes, n_edges, edges) = io::graph(io::lines(opt.instance()));

    e =  BoolVarArray(*this, (int)edges.size(), 0, 1);

    spanning_tree_ctr(*this, e, edges, n_nodes);

//    cbsbranch(*this, e,  CBSBranchingHeuristic::MAX_SD);
    branch(*this, e, BOOL_VAR_NONE(), BOOL_VAL_MIN());
  }

  SpanningTree(bool share, SpanningTree& s)
    : Script(share,s) {
    e.update(*this, share, s.e);
  }

  Space* copy(bool share) override {
    return new SpanningTree(share,*this);
  }

  void print(std::ostream& os) const override {
    for (int i=0; i<e.size(); i++) {
     os << e[i] << std::endl;
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
