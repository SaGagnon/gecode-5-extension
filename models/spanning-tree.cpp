#include <set>
#include <stack>
#include <armadillo>

#include <gecode/driver.hh>
#include <unordered_map>

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
using NaryProp = NaryPropagator<BoolView, Int::ME_INT_DOM>;

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

    // unvisited nodes
    auto u_nodes = utils::set_zero_to<node_t>(graph.n_nodes);
    while (!u_nodes.empty()) {
      std::set<node_t> CC;

      // Construction of CC (connected component) without cycles
      auto ret = DFS_one(*u_nodes.begin(), [&](node_t node) {
        if (!CC.insert(node).second)
          return ES_FAILED;
        u_nodes.erase(node);
        return ES_OK;
      });
      if (ret != ES_OK) return ret;

      // If there's only one node
      if (CC.size() == 1) {

        // Number of adjacent unassigned edges
        int n_unassigned = 0;
        BoolView last_view_unassigned;

        ret = DFS_one(*CC.begin(), [&](node_t node) {
          for_each_adj(node, [&](node_t adj_n, BoolView& adj_v ) {
            if (adj_v.none()) {
              n_unassigned++;
              last_view_unassigned = adj_v;
            }
          });
          return ES_OK;
        });
        if (ret != ES_OK) return ret;

        // If there's no adjacent edges, we fail
        if (n_unassigned == 0)
          return ES_FAILED;
        // If there's only one, we are obliged to take it
        if (n_unassigned == 1)
          last_view_unassigned.eq(home, 1);

      } else if (CC.size() > 3) {

        // We make sure there's no unassigned edges pointing to the CC itself
        ret = DFS_one(*CC.begin(), [&](node_t node) {
          for_each_adj(node, [&](node_t adj_n, BoolView& adj_v ) {
            if (adj_v.none() && CC.find(adj_n) != CC.end())
              adj_v.eq(home, 0);
          });
          return ES_OK;
        });
        if (ret != ES_OK) return ret;
      }

      // If the connected component does not contain the whole graph, it must
      // at least include an unassigned edge (for the whole graph to be
      // connected)
      if (CC.size() != graph.n_nodes) {
        unsigned int n_unassigned = 0;
        ret = DFS_one(*CC.begin(),
                   [&](node_t node) {
                     for_each_adj(node, [&](node_t adj_n, BoolView& adj_v ) {
                       if (adj_v.none())
                         n_unassigned += 1;
                     });
                     return ES_OK;
                   });
        if (ret != ES_OK) return ret;

        if (n_unassigned == 0)
          return ES_FAILED;
      }
    }

    return ES_FIX;
  }

  void solndistrib(Space& home, SolnDistrib *dist) const override {

    // each node belong to a connected component identified by a master node
    std::unordered_map<node_t, node_t> ncc;
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
      laplacian = arma::mat(graph.n_nodes, graph.n_nodes, arma::fill::zeros);
      DFS_all(0, [&](node_t node) {
        if (ncc[node] != node)
          laplacian(node, node) = 1;
        for_each_adj(node, [&](node_t adj_n, BoolView& adj_v) {
          // We are going to see each edge two times and we only want
          // to add it one time
          if (adj_v.none() && node < adj_n) {
            auto cc1 = ncc[node];
            auto cc2 = ncc[adj_n];
            assert(cc1 != cc2);
            laplacian(cc1, cc2) -=1;
            laplacian(cc2, cc1) -=1;
            laplacian(cc1, cc1) += 1;
            laplacian(cc2, cc2) += 1;
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
//      rel(*this, sum(adj_edges[i]) <= 2);
    }

    cbsbranch(*this, e,  CBSBranchingHeuristic::MAX_SD);
//    branch(*this, e, BOOL_VAR_NONE(), BOOL_VAL_MIN());
  }

  ConstrainedSpanningTree(bool share, ConstrainedSpanningTree& s)
    : Script(share,s), instance(s.instance) {
    e.update(*this, share, s.e);
  }

  Space* copy(bool share) override {
    return new ConstrainedSpanningTree(share,*this);
  }

  void print(std::ostream& os) const override {
    bool finished = true;
    for (int i=0; i<e.size(); i++)
      if (e[i].none()) finished = false;

    if (finished) {
      n_nodes_t           n_nodes;
      n_edges_t           n_edges;
      std::vector<edge_t> edges;

      std::tie(n_nodes, n_edges, edges) = io::graph(io::lines(instance));

      std::vector<std::vector<node_t>> adj_list(n_nodes);
      {
        n_edges_t n_edges_sol = 0;
        for (int i = 0; i < n_edges; i++) {
          if (e[i].one()) {
            node_t n1, n2;
            std::tie(n1, n2) = edges[i];
            adj_list[n1].push_back(n2);
            adj_list[n2].push_back(n1);

            n_edges_sol++;
          }
        }
        assert(n_edges_sol == n_nodes - 1);
      }


      for (int i=0; i<n_nodes; i++)
        assert(!adj_list[i].empty());


      auto u_nodes = utils::set_zero_to<node_t>(n_nodes);

      node_t start = *u_nodes.begin();

      using parent_node_t = node_t;
      std::stack<std::pair<parent_node_t, node_t>> stack;
      stack.emplace(start, start);

      while (!stack.empty()) {
        parent_node_t parent; node_t node;
        std::tie(parent, node) = stack.top();
        stack.pop();

        auto n_ereased = u_nodes.erase(node);
        assert(n_ereased == 1);

        for (unsigned int adj_n : adj_list[node]) {
          if (adj_n != parent)
            stack.emplace(node, adj_n);
        }
      }
      assert(u_nodes.empty());
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

  Script::run<ConstrainedSpanningTree,DFS,InstanceOptions>(opt);
}
