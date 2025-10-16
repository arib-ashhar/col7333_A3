#include <bits/stdc++.h>
using namespace std;

struct Point { int x{}, y{}; };
struct MetroLineIO { Point s, e; };

struct Scenario1IO {
    int N{}, M{}, K{}, J{};
    vector<MetroLineIO> lines;
};

struct Scenario2IO {
    int N{}, M{}, K{}, J{}, P{};
    vector<MetroLineIO> lines;
    vector<Point> popular;
};

static bool in_bounds(int x, int y, int N, int M) {
    return (0 <= x && x < N && 0 <= y && y < M);
}

static inline int manhattan(int x1,int y1,int x2,int y2){
    return abs(x1-x2)+abs(y1-y2);
}

static Scenario1IO parseScenario1(istream& in) {
    Scenario1IO sc;
    if (!(in >> sc.N >> sc.M >> sc.K >> sc.J)) {
        throw runtime_error("Failed to read N M K J for Scenario 1.");
    }
    if (sc.N <= 0 || sc.M <= 0) {
        throw runtime_error("N and M must be positive.");
    }
    if (sc.K < 0 || sc.J < 0) {
        throw runtime_error("K and J must be non-negative.");
    }

    sc.lines.resize(sc.K);
    // Checking uniqueness of start/end locations
    auto key = [](int x, int y){ return ( (uint64_t(x) << 32) ^ uint32_t(y) ); };
    unordered_set<uint64_t> endpoints;
    endpoints.reserve(size_t(sc.K) * 2u);

    for (int k = 0; k < sc.K; ++k) {
        int sx, sy, ex, ey;
        if (!(in >> sx >> sy >> ex >> ey)) {
            throw runtime_error("Failed to read (sx sy ex ey) for line " + to_string(k));
        }
        if (!in_bounds(sx, sy, sc.N, sc.M) || !in_bounds(ex, ey, sc.N, sc.M)) {
            ostringstream oss;
            oss << "Line " << k << " has out-of-bounds endpoint(s): "
                << "(" << sx << "," << sy << ") -> (" << ex << "," << ey << ") "
                << "but grid is 0.." << sc.N-1 << " x 0.." << sc.M-1;
            throw runtime_error(oss.str());
        }
        uint64_t ks = key(sx, sy), ke = key(ex, ey);
        if (!endpoints.insert(ks).second) {
            ostringstream oss; oss << "Duplicate endpoint at (" << sx << "," << sy << ")";
            throw runtime_error(oss.str());
        }
        if (!endpoints.insert(ke).second) {
            ostringstream oss; oss << "Duplicate endpoint at (" << ex << "," << ey << ")";
            throw runtime_error(oss.str());
        }
        sc.lines[k] = MetroLineIO{ Point{sx, sy}, Point{ex, ey} };
    }
    return sc;
}

static Scenario2IO parseScenario2(istream& in) {
    Scenario2IO sc;
    if (!(in >> sc.N >> sc.M >> sc.K >> sc.J >> sc.P))
        throw runtime_error("Failed to read N M K J P for Scenario 2.");
    if (sc.N <= 0 || sc.M <= 0) throw runtime_error("N and M must be positive.");
    if (sc.K < 0 || sc.J < 0 || sc.P < 0) throw runtime_error("K, J, P must be non-negative.");

    sc.lines.resize(sc.K);
    auto key = [](int x, int y){ return ( (uint64_t(x) << 32) ^ uint32_t(y) ); };
    unordered_set<uint64_t> endpoints; endpoints.reserve(size_t(sc.K) * 2u);

    for (int k = 0; k < sc.K; ++k) {
        int sx, sy, ex, ey;
        if (!(in >> sx >> sy >> ex >> ey))
            throw runtime_error("Failed to read (sx sy ex ey) for line " + to_string(k));
        if (!in_bounds(sx, sy, sc.N, sc.M) || !in_bounds(ex, ey, sc.N, sc.M))
            throw runtime_error("Scenario 2: endpoint out-of-bounds");
        uint64_t ks = key(sx, sy), ke = key(ex, ey);
        if (!endpoints.insert(ks).second || !endpoints.insert(ke).second)
            throw runtime_error("Scenario 2: duplicate endpoint");
        sc.lines[k] = MetroLineIO{ Point{sx, sy}, Point{ex, ey} };
    }

    sc.popular.resize(sc.P);
    for (int p = 0; p < sc.P; ++p) {
        int x, y;
        if (!(in >> x >> y)) throw runtime_error("Failed to read popular cell");
        if (!in_bounds(x, y, sc.N, sc.M)) throw runtime_error("Popular cell out-of-bounds");
        sc.popular[p] = Point{x, y};
    }
    return sc;
}

//------------------CNF + Encoding ------------------------
struct CNF {
    int maxVar = 0;
    vector<vector<int>> clauses;

    int newVar() { return ++maxVar; }

    void addClause(const initializer_list<int>& lits) { clauses.emplace_back(lits); }
    void addClause(const vector<int>& lits)          { clauses.push_back(lits); }

    // Helper for Boolean expressions into CNF form
    // a -> b = (-a, b)
    void imply(int a, int b) { addClause({-a, b}); }

    // z <-> (a AND b) = (-z, a) ^ (-z, b) ^ (-a, -b, z)
    void eq_and(int z, int a, int b) {
        addClause({-z, a});
        addClause({-z, b});
        addClause({-a, -b, z});
    }

    // z <-> (OR S)
    int eq_or(const vector<int>& S) {
        int z = newVar();
        if (S.empty()) {
            // z = false
            addClause({-z});
            return z;
        }
        // (¬s_i ∨ z) for all i
        for (int s : S) addClause({-s, z});
        // (¬z ∨ s1 ∨ s2 ∨ ... )
        {
            vector<int> c;
            c.reserve(S.size()+1);
            c.push_back(-z);
            c.insert(c.end(), S.begin(), S.end());
            addClause(c);
        }
        return z;
    }

    // AtLeastOne(S)
    void atLeastOne(const vector<int>& S) {
        if (!S.empty()) addClause(S);
    }

    // Pairwise AtMostOne(S)
    void atMostOne_pairwise(const vector<int>& S) {
        for (size_t i = 0; i < S.size(); ++i)
            for (size_t j = i+1; j < S.size(); ++j)
                addClause({-S[i], -S[j]});
    }

    // ExactlyOne(S) = AtLeastOne + AtMostOne
    void exactlyOne(const vector<int>& S) {
        atLeastOne(S);
        atMostOne_pairwise(S);
    }

    // exactlyOne is true when S1 is true or else all false: (S1 -> ExactlyOne(S))
    void exactlyOne_GivenS1(int S1, const vector<int>& S) {
        if (S.empty()) return;
        {
            vector<int> c; 
            c.reserve(S.size()+1);

            // If S1 is true, exactly one statement in S is true
            c.push_back(-S1);
            c.insert(c.end(), S.begin(), S.end());
            addClause(c); // (-S1 ∨ OR S)
        }
       // If S1 is false, all other pair  (i,j) in S becomes irrrelevant
        for (size_t i = 0; i < S.size(); ++i)
            for (size_t j = i+1; j < S.size(); ++j)
                addClause({-S1, -S[i], -S[j]});
    }

    // Guarded AllZero: (¬S1 -> all s_i = 0)  == (S1 ∨ ¬s_i)
    void allZero_Given_S1_False(const vector<int>& S, int S1) {
        for (int s : S) addClause({S1, -s});
    }

    // Cardinality clause, will be used to add clause for J bends in a metro Line
    // Enforce sum(P) <= K using Sinz sequential counter.
    void atMostK(const vector<int>& P, int K) {
        const int n = (int)P.size();
        if (K < 0) return;
        if (K == 0) { for (int p : P) addClause({-p}); return; }
        if (K >= n) return;

        unordered_map<long long,int> memo;
        memo.reserve((n-1)*max(1,K));

        auto s = [&](int i, int j)->int {
            long long key = ((long long)i<<20) ^ j ^ 0x9e3779b97f4a7c15ULL;
            auto it = memo.find(key);
            if (it != memo.end()) return it->second;
            int v = newVar();
            memo.emplace(key, v);
            return v;
        };

        addClause({-P[0], s(1,1)});
        for (int j = 2; j <= K; ++j) 
            addClause({-s(1,j)});
        for (int i = 2; i <= n-1; ++i) {
            addClause({-P[i-1], s(i,1)});
            addClause({-s(i-1,1), s(i,1)});
            for (int j = 2; j <= K; ++j) {
                addClause({-s(i-1,j), s(i,j)});
                addClause({-P[i-1], -s(i-1,j-1), s(i,j)});
            }
        }
        for (int i = 2; i <= n; ++i) {
            addClause({-P[i-1], -s(i-1, K)});
        }
    }

    void writeDimacsIntoSATinput(const string& path) {
        ofstream out(path, ios::out|ios::trunc);
        if (!out) throw runtime_error("Could not open " + path + " for writing.");
        out << "p cnf " << maxVar << " " << clauses.size() << "\n";
        for (auto& c : clauses) {
            for (int lit : c) out << lit << " ";
            out << "0\n";
        }
    }
};

// ---------------- Grid indexing & edges ----------------------
struct Grid {
    int N, M; // width (x), height (y)
    // node id: u = y*N + x
    int node(int x, int y) const { return y*N + x; }

    // Directed edges list e: (u->v)
    vector<pair<int,int>> edges;
    vector<vector<int>> out_adj; // u -> [edge indices]
    vector<vector<int>> in_adj;  // v -> [edge indices]
    vector<bool> isHoriz;        // edges[e] is horizontal

    Grid(int N_, int M_): N(N_), M(M_) {
        int V = N*M;
        out_adj.assign(V, {});
        in_adj.assign(V, {});
        // Build edges (4-neighbour)
        for (int y=0; y<M; ++y) 
            for (int x=0; x<N; ++x) {
                int u = node(x,y);
                const int dx[4]={1,-1,0,0}, dy[4]={0,0,1,-1};
                for (int d=0; d<4; ++d) {
                    int nx=x+dx[d], ny=y+dy[d];
                    if (!in_bounds(nx, ny, N, M)) continue;
                    int v = node(nx,ny);
                    int idx = (int)edges.size();
                    edges.push_back({u,v});
                    out_adj[u].push_back(idx);
                    in_adj[v].push_back(idx);
                }
            }
        isHoriz.resize(edges.size());
        for (size_t i=0;i<edges.size();++i) {
            auto [u,v]=edges[i];
            int uy=u/N, vy=v/N;
            isHoriz[i] = (uy==vy);
        }
    }
};

// ------------- Variable manager (X, F, D, T) ----------------
struct VarMgr {
    CNF& cnf;
    ofstream* vmap;
    int N, M, Umax, K; // Umax = max level per line (inclusive)
    VarMgr(CNF& c, int N_, int M_, int K_, int Umax_)
        : cnf(c), N(N_), M(M_), Umax(Umax_), K(K_) {}

    static uint64_t pack4(int a,int b,int c,int d){
        return ( (uint64_t)(uint32_t)a << 48 )
             ^ ( (uint64_t)(uint32_t)b << 32 )
             ^ ( (uint64_t)(uint32_t)c << 16 )
             ^ ( (uint64_t)(uint32_t)d );
    }
    static uint64_t pack3(int a,int b,int c){
        return ( (uint64_t)(uint32_t)a << 40 )
             ^ ( (uint64_t)(uint32_t)b << 20 )
             ^ ( (uint64_t)(uint32_t)c );
    }

    unordered_map<uint64_t,int> Xvar;     // X(k,u)  - Line k uses cell u on the grid
    unordered_map<uint64_t,int> Fvar;     // F(k, edgeIndex)  - Line k used edge with id=edgeIndex
    unordered_map<uint64_t,int> Dvar;     // D(k,u,t)   -  cell u is at distance t for Line k, where t = [1....U-1]
    unordered_map<uint64_t,int> Tvar;     // T(k,u)   - Line k turns at cell u

    int var_X(int k, int u) {
        uint64_t key = pack3(1,k,u);
        auto it = Xvar.find(key);
        if (it!=Xvar.end()) return it->second;
        int v = cnf.newVar();
        Xvar[key]=v; 
        if (vmap) (*vmap) << "X " << k << " " << u << " " << v << "\n";
        return v;
    }

    int var_F(int k, int eidx) {
        uint64_t key = pack3(2,k,eidx);
        auto it = Fvar.find(key);
        if (it!=Fvar.end()) return it->second;
        int v = cnf.newVar();
        Fvar[key]=v; 
        if (vmap) (*vmap) << "F " << k << " " << eidx << " " << v << "\n";
        return v;
    }

    int var_D(int k, int u, int t) {
        uint64_t key = pack4(3,k,u,t);
        auto it = Dvar.find(key);
        if (it!=Dvar.end()) return it->second;
        int v = cnf.newVar();
        Dvar[key]=v; 
        if (vmap) (*vmap) << "D " << k << " " << u << " " << t << " " << v << "\n";
        return v;
    }

    int var_T(int k, int u) {
        uint64_t key = pack3(4,k,u);
        auto it = Tvar.find(key);
        if (it!=Tvar.end()) return it->second;
        int v = cnf.newVar();
        Tvar[key]=v; 
        if (vmap) (*vmap) << "T " << k << " " << u << " " << v << "\n";
        return v;
    }
};

// -------------------- Encoders per line ----------------------

// Edge implies endpoints used; also forbid two-way on same undirected edge
static void encode_edge_endpoint_and_no_2way(int k, const Grid& G, VarMgr& vm) {
    // Edge -> endpoints used
    for (int e = 0; e < (int)G.edges.size(); ++e) {
        auto [u,v] = G.edges[e];
        int F = vm.var_F(k, e);
        int Xu = vm.var_X(k, u);
        int Xv = vm.var_X(k, v);
        vm.cnf.imply(F, Xu);
        vm.cnf.imply(F, Xv);
    }
    // No 2-way on same undirected edge
    unordered_map<long long,int> rev;
    rev.reserve(G.edges.size()*2);
    auto key = [&](int a,int b)->long long { return ( (long long)a<<32 ) ^ (long long)b; };
    for (int e=0;e<(int)G.edges.size();++e) {
        int u=G.edges[e].first, v=G.edges[e].second;
        rev[key(v,u)] = e;
    }
    vector<char> done(G.edges.size(), 0);
    for (int e=0;e<(int)G.edges.size();++e) if (!done[e]) {
        int u=G.edges[e].first, v=G.edges[e].second;
        auto it=rev.find(key(u,v));
        if (it!=rev.end()) {
            int r = it->second;
            if (r==e) { done[e]=1; continue; }
            int Fe = vm.var_F(k, e);
            int Fr = vm.var_F(k, r);
            vm.cnf.addClause({-Fe, -Fr});
            done[e]=done[r]=1;
        }
    }
}

// Degree constraints at nodes:
// start s: exactlyOne out, all incoming = 0, X true
// end   e: exactlyOne in,  all outgoing = 0, X true
// others u: (X -> exactlyOne in) & (X -> exactlyOne out) & ((¬X) -> all incident F=0)
static void encode_degrees(int k, const Grid& G, VarMgr& vm, int s_u, int e_u) {
    for (int u=0; u<(int)G.out_adj.size(); ++u) {
        int X = vm.var_X(k, u);
        const auto& OUT = G.out_adj[u];
        const auto& IN  = G.in_adj[u];
        vector<int> Fout, Fin;
        Fout.reserve(OUT.size());
        Fin.reserve(IN.size());
        for (int e : OUT) Fout.push_back(vm.var_F(k,e));
        for (int e : IN)  Fin.push_back(vm.var_F(k,e));

        if (u == s_u) {
            vm.cnf.exactlyOne(Fout);
            for (int f : Fin) vm.cnf.addClause({-f});
            vm.cnf.addClause({X});
        } else if (u == e_u) {
            vm.cnf.exactlyOne(Fin);
            for (int f : Fout) vm.cnf.addClause({-f});
            vm.cnf.addClause({X});
        } else {
            vm.cnf.exactlyOne_GivenS1(X, Fin);
            vm.cnf.exactlyOne_GivenS1(X, Fout);
            vector<int> all = Fin; all.insert(all.end(), Fout.begin(), Fout.end());
            vm.cnf.allZero_Given_S1_False(all, X);
        }
    }
}

// Level constraints:
// intially, D(s,0)=true; D(u!=s,0)=false
// X(u) -> ExactlyOne_t D(u,t) ; (¬X(u)) -> all D(u,t)=0
// F(u->v) ∧ D(u,t) -> D(v,t+1)
// D(u,t) -> OR_w (F(w->u) ∧ D(w,t-1))
static void encode_levels(int k, const Grid& G, VarMgr& vm, int s_u, int e_u, int U) {
    CNF& cnf = vm.cnf;

    // Anchor t=0
    cnf.addClause({ vm.var_D(k, s_u, 0) });
    for (int u=0; u<(int)G.out_adj.size(); ++u) if (u != s_u)
        cnf.addClause({ -vm.var_D(k, u, 0) });

    // X gating of D
    for (int u=0; u<(int)G.out_adj.size(); ++u) {
        int X = vm.var_X(k,u);
        vector<int> Dt; Dt.reserve(U+1);
        for (int t=0; t<=U; ++t) Dt.push_back(vm.var_D(k,u,t));
        cnf.exactlyOne_GivenS1(X, Dt);
        for (int t=0; t<=U; ++t) cnf.addClause({ X, -Dt[t] });
    }

    // F(u->v) ∧ D(u,t) -> D(v,t+1)
    for (int e=0; e<(int)G.edges.size(); ++e) {
        auto [u,v]=G.edges[e];
        int F = vm.var_F(k,e);
        for (int t=0; t<U; ++t) {
            int Du = vm.var_D(k,u,t);
            int Dv = vm.var_D(k,v,t+1);
            cnf.addClause({ -F, -Du, Dv });
        }
    }

    // Predecessor existence
    for (int u=0; u<(int)G.in_adj.size(); ++u) {
        if (u == s_u) continue;
        const auto& IN = G.in_adj[u];
        for (int t=1; t<=U; ++t) {
            int Du_t = vm.var_D(k,u,t);
            vector<int> Ps; Ps.reserve(IN.size());
            for (int e_in : IN) {
                int w   = G.edges[e_in].first;
                int Fwu = vm.var_F(k, e_in);
                int Dw  = vm.var_D(k, w, t-1);
                int P   = cnf.newVar();
                cnf.eq_and(P, Fwu, Dw);
                Ps.push_back(P);
            }
            vector<int> clause; clause.reserve(Ps.size()+1);
            clause.push_back(-Du_t);
            clause.insert(clause.end(), Ps.begin(), Ps.end());
            cnf.addClause(clause);
        }
    }
}

// Turn indicator T(u) and AtMost(J) over all u != s,e
// T ↔ ((HIn ∧ VOut) ∨ (VIn ∧ HOut)), and T -> X
static void encode_turns_and_limitJ(int k, const Grid& G, VarMgr& vm, int s_u, int e_u, int J) {
    CNF& cnf = vm.cnf;
    vector<int> Tvars;

    for (int u=0; u<(int)G.out_adj.size(); ++u) {
        int T = vm.var_T(k,u);
        int X = vm.var_X(k,u);

        if (u==s_u || u==e_u) {
            cnf.addClause({-T}); // endpoints not turns
            continue;
        }

        // Collect incoming/outgoing edges by orientation
        vector<int> HIn, VIn, HOut, VOut;
        for (int e : G.in_adj[u])  (G.isHoriz[e] ? HIn : VIn).push_back(vm.var_F(k,e));
        for (int e : G.out_adj[u]) (G.isHoriz[e] ? HOut: VOut).push_back(vm.var_F(k,e));

        int hIn  = HIn.empty()  ? 0 : cnf.eq_or(HIn);
        int vIn  = VIn.empty()  ? 0 : cnf.eq_or(VIn);
        int hOut = HOut.empty() ? 0 : cnf.eq_or(HOut);
        int vOut = VOut.empty() ? 0 : cnf.eq_or(VOut);

        auto mk_and = [&](int a,int b)->int {
            if (a==0 || b==0) { int z=cnf.newVar(); cnf.addClause({-z}); return z; }
            int z=cnf.newVar(); cnf.eq_and(z,a,b); return z;
        };
        int A = mk_and(hIn, vOut);
        int B = mk_and(vIn, hOut);

        int OrAB = cnf.eq_or({A,B});
        // T ↔ OrAB
        cnf.addClause({-T, OrAB});
        cnf.addClause({-OrAB, T});

        // T -> X
        cnf.imply(T, X);

        Tvars.push_back(T);
    }

    cnf.atMostK(Tvars, J);
}

// Cell exclusivity across lines (at most one line per location)
static void encode_cell_exclusivity_allK(const Grid& G, VarMgr& vm, int K, const unordered_set<int>& popular_set) {
    int V = G.N * G.M;
    for (int u = 0; u < V; ++u) {
        if (popular_set.count(u)) continue; // allow multiple lines on popular cells
        for (int k = 0; k < K; ++k) {
            int Xu = vm.var_X(k, u);
            for (int k2 = k+1; k2 < K; ++k2) {
                int Xu2 = vm.var_X(k2, u);
                vm.cnf.addClause({-Xu, -Xu2});
            }
        }
    }
}

// Allow multiple metro lines to pass through one popular cell
static void encode_popular_cell_for_allK(const vector<int>& popular_nodes, int K, VarMgr& vm, CNF& cnf){
    for (int u_pop : popular_nodes) {
        vector<int> orX; orX.reserve(K);
        for (int k = 0; k < K; ++k) 
            orX.push_back(vm.var_X(k, u_pop));
        cnf.addClause(orX);
    }
}

struct MetroLine { int sx, sy, ex, ey; }; // compact carrier
struct Scenario1Data {
    int N, M, K, J;
    vector<MetroLine> lines;
};

static void encodeScenario1(const Scenario1Data& sc, const string& outPath, const string& vrmapPath) {
    Grid G(sc.N, sc.M);
    CNF cnf;
    ofstream vmap(vrmapPath, ios::out | ios::trunc);
    if (!vmap) throw runtime_error("Cannot open " + vrmapPath);

    // Conservative level bound (inclusive): 0..N*M-1
    int Umax = sc.N*sc.M - 1;
    VarMgr vm(cnf, sc.N, sc.M, sc.K, Umax);
    vm.vmap = &vmap;

    // Per-line constraints
    for (int k=0; k<sc.K; ++k) {
        int s = G.node(sc.lines[k].sx, sc.lines[k].sy);
        int e = G.node(sc.lines[k].ex, sc.lines[k].ey);
        int U = min(sc.N*sc.M - 1, manhattan(sc.lines[k].sx, sc.lines[k].sy, sc.lines[k].ex, sc.lines[k].ey) + 3);

        encode_edge_endpoint_and_no_2way(k, G, vm);
        encode_degrees(k, G, vm, s, e);
        encode_levels(k, G, vm, s, e, U);
        encode_turns_and_limitJ(k, G, vm, s, e, sc.J);
    }

    // Inter-line exclusivity
    encode_cell_exclusivity_allK(G, vm, sc.K, {});
    vmap.close();

    // Write DIMACS
    cnf.writeDimacsIntoSATinput(outPath);
}

static void encodeScenario2(const Scenario1Data& sc, const vector<Point>& popular_points, const string& outPath, const string& vrmapPath) {
    Grid G(sc.N, sc.M);
    CNF cnf;
    ofstream vmap(vrmapPath, ios::out | ios::trunc);
    if (!vmap) throw runtime_error("Cannot open " + vrmapPath);

    int Umax = sc.N*sc.M - 1;
    VarMgr vm(cnf, sc.N, sc.M, sc.K, Umax);
    vm.vmap = &vmap;

    // Precompute set of popular node-ids
    vector<int> popular_nodes; popular_nodes.reserve(popular_points.size());
    unordered_set<int> popular_set;
    for (auto &p : popular_points) {
        int u = G.node(p.x, p.y);
        popular_nodes.push_back(u);
        popular_set.insert(u);
    }

    // Per-line constraints
    for (int k = 0; k < sc.K; ++k) {
        int s = G.node(sc.lines[k].sx, sc.lines[k].sy);
        int e = G.node(sc.lines[k].ex, sc.lines[k].ey);
        int U = min(sc.N*sc.M - 1, manhattan(sc.lines[k].sx, sc.lines[k].sy, sc.lines[k].ex, sc.lines[k].ey) + 3);
        encode_edge_endpoint_and_no_2way(k, G, vm);
        encode_degrees(k, G, vm, s, e);
        encode_levels(k, G, vm, s, e, U);
        encode_turns_and_limitJ(k, G, vm, s, e, sc.J);
    }

    encode_popular_cell_for_allK(popular_nodes, sc.K, vm, cnf);
    encode_cell_exclusivity_allK(G, vm, sc.K, popular_set);

    vmap.close();
    cnf.writeDimacsIntoSATinput(outPath);
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <base>\n"
             << "Example: " << argv[0] << " foo   # reads foo.city, writes foo.satinput\n";
        return 1;
    }
    const string base = argv[1];
    const string inPath  = base + ".city";
    const string outPath = base + ".satinput";
    const string vrmapPath = base + ".varmap";

    ifstream fin(inPath);
    if (!fin) {
        cerr << "Error: cannot open file '" << inPath << "'\n";
        return 1;
    }

    int scenarioType = 0;
    if (!(fin >> scenarioType)) {
        cerr << "Error: failed to read scenario type (expected 1 or 2).\n";
        return 1;
    }

    try {
        switch (scenarioType) {
            case 1: {
                Scenario1IO scio = parseScenario1(fin);

                Scenario1Data data;
                data.N = scio.N;
                data.M = scio.M;
                data.K = scio.K;
                data.J = scio.J;
                data.lines.reserve(scio.K);
                for (int k = 0; k < scio.K; ++k) {
                    data.lines.push_back({
                        scio.lines[k].s.x, scio.lines[k].s.y,
                        scio.lines[k].e.x, scio.lines[k].e.y
                    });
                }

                encodeScenario1(data, outPath, vrmapPath);

                cout << "Wrote DIMACS CNF to " << outPath << "\n";
                cout << "Grid " << data.N << "x" << data.M
                     << ", K=" << data.K << ", J=" << data.J << "\n";
                break;
            }
            case 2: {
                Scenario2IO sc2 = parseScenario2(fin);

                Scenario1Data data;
                data.N = sc2.N;
                data.M = sc2.M;
                data.K = sc2.K;
                data.J = sc2.J;
                data.lines.reserve(sc2.K);
                for (int k = 0; k < sc2.K; ++k) {
                    data.lines.push_back({
                        sc2.lines[k].s.x, sc2.lines[k].s.y,
                        sc2.lines[k].e.x, sc2.lines[k].e.y
                    });
                }

                vector<Point> popular = sc2.popular;
                encodeScenario2(data, popular, outPath, vrmapPath);

                cout << "Scenario 2: wrote DIMACS CNF to " << outPath << "\n";
                cout << "Grid " << data.N << "x" << data.M
                     << ", K=" << data.K << ", J=" << data.J << ", P=" << sc2.P << "\n";
                break;
            }
            default:
                cerr << "Error: unsupported scenario type '" << scenarioType
                     << "'. Expected 1 or 2.\n";
                return 1;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}