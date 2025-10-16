#include <bits/stdc++.h>
using namespace std;

// ---------- Basic IO structs (same as encoder) ----------
struct Point { int x{}, y{}; };
struct MetroLine { Point s, e; };

struct Scenario1 {
    int N{}, M{}, K{}, J{};
    vector<MetroLine> lines;
};

struct Scenario2 {
    int N{}, M{}, K{}, J{}, P{};
    vector<MetroLine> lines;
    vector<Point> popular;
};

static bool in_bounds(int x, int y, int N, int M) {
    return (0 <= x && x < N && 0 <= y && y < M);
}

static Scenario1 parseScenario1(istream& in) {
    Scenario1 sc;
    if (!(in >> sc.N >> sc.M >> sc.K >> sc.J)) {
        throw runtime_error("Failed to read N M K J for Scenario 1.");
    }
    if (sc.N <= 0 || sc.M <= 0) throw runtime_error("N and M must be positive.");
    if (sc.K < 0 || sc.J < 0) throw runtime_error("K and J must be non-negative.");

    sc.lines.resize(sc.K);
    auto key = [](int x, int y){ return ( (uint64_t(x) << 32) ^ uint32_t(y) ); };
    unordered_set<uint64_t> endpoints;
    endpoints.reserve(size_t(sc.K) * 2u);

    for (int k = 0; k < sc.K; ++k) {
        int sx, sy, ex, ey;
        if (!(in >> sx >> sy >> ex >> ey))
            throw runtime_error("Failed to read (sx sy ex ey) for line " + to_string(k));
        if (!in_bounds(sx, sy, sc.N, sc.M) || !in_bounds(ex, ey, sc.N, sc.M)) {
            ostringstream oss;
            oss << "Line " << k << " endpoints out of bounds: "
                << "(" << sx << "," << sy << ") -> (" << ex << "," << ey << ")";
            throw runtime_error(oss.str());
        }
        uint64_t ks = key(sx, sy), ke = key(ex, ey);
        if (!endpoints.insert(ks).second) throw runtime_error("Duplicate start endpoint");
        if (!endpoints.insert(ke).second) throw runtime_error("Duplicate end endpoint");
        sc.lines[k] = MetroLine{ Point{sx, sy}, Point{ex, ey} };
    }
    return sc;
}

static Scenario2 parseScenario2(istream& in) {
    Scenario2 sc;
    if (!(in >> sc.N >> sc.M >> sc.K >> sc.J >> sc.P)) {
        throw runtime_error("Failed to read N M K J P for Scenario 2.");
    }
    if (sc.N <= 0 || sc.M <= 0) throw runtime_error("N and M must be positive.");
    if (sc.K < 0 || sc.J < 0 || sc.P < 0) throw runtime_error("K, J, P must be non-negative.");

    sc.lines.resize(sc.K);
    for (int k = 0; k < sc.K; ++k) {
        int sx, sy, ex, ey;
        if (!(in >> sx >> sy >> ex >> ey))
            throw runtime_error("Failed to read (sx sy ex ey) for line " + to_string(k));
        if (!in_bounds(sx, sy, sc.N, sc.M) || !in_bounds(ex, ey, sc.N, sc.M))
            throw runtime_error("Scenario 2: endpoint out-of-bounds");
        sc.lines[k] = MetroLine{ Point{sx, sy}, Point{ex, ey} };
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

// ---------- Grid to translate node ids based on (x,y) ----------
struct Grid {
    int N, M;
    Grid(int N_, int M_) : N(N_), M(M_) {}
    inline int node(int x, int y) const { return y*N + x; }
    inline pair<int,int> xy(int u) const { return {u%N, u/N}; }
};

// ---------- MiniSAT model reader ----------
static bool readModel(const string& path, vector<char>& value, int& maxVar) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open model file: " + path);
    string token;
    maxVar = 0;
    bool sat = false;

    vector<int> lits;
    while (in >> token) {
        if (token == "UNSAT") { sat = false; break; }
        if (token == "SAT")   { sat = true;  continue; }
        if (token == "v")     { 
            int lit;
            while (in >> lit) {
                if (lit == 0) break;
                lits.push_back(lit);
                maxVar = max(maxVar, abs(lit));
            }
        } else {
            int lit;
            try { lit = stoi(token); }
            catch(...) { continue; }
            if (lit == 0) continue;
            lits.push_back(lit);
            maxVar = max(maxVar, abs(lit));
        }
    }

    if (!sat) return false;
    value.assign(maxVar + 1, 0);
    for (int lit : lits) {
        int v = abs(lit);
        value[v] = (lit > 0) ? 1 : -1;
    }
    return true;
}

// ---------- Varmap reader ----------
// Format (whitespace-separated, one per line):
//   D k u t var
//   F k e var
//   X k u var
//   T k u var
struct VarMap {
    // We only need D to reconstruct the path. (F/X/T are optional here.)
    struct DKey { int k,u,t; };
    struct DKeyHash {
        size_t operator()(DKey const& a) const noexcept {
            return ( (size_t)a.k * 1315423911u ) ^ ((size_t)a.u<<1) ^ (size_t)a.t;
        }
    };
    struct DKeyEq {
        bool operator()(DKey const& a, DKey const& b) const noexcept {
            return a.k==b.k && a.u==b.u && a.t==b.t;
        }
    };
    unordered_map<DKey,int,DKeyHash,DKeyEq> Dvar;

    void load(const string& path) {
        ifstream in(path);
        if (!in) throw runtime_error("Cannot open varmap file: " + path);
        string typ;
        while (in >> typ) {
            if (typ == "D") {
                int k,u,t,var; in >> k >> u >> t >> var;
                Dvar[{k,u,t}] = var;
            } else if (typ == "F") {
                int k,e,var; in >> k >> e >> var; (void)k; (void)e; (void)var;
                // Not needed for decoding directions
            } else if (typ == "X") {
                int k,u,var; in >> k >> u >> var; (void)k; (void)u; (void)var;
            } else if (typ == "T") {
                int k,u,var; in >> k >> u >> var; (void)k; (void)u; (void)var;
            } else {
                // skip line
                string rest; getline(in, rest);
            }
        }
    }
};

// ---------- Reconstruct per-line path from D-variables ----------
static vector<vector<char>> reconstructPaths(
    const Scenario1& sc,
    const Grid& G,
    const VarMap& vm,
    const vector<char>& val // index by var id: 1 -> true, -1 -> false, 0 -> unknown
){
    vector<vector<char>> paths(sc.K);

    for (int k = 0; k < sc.K; ++k) {
        // Find all (t,u) such that D_{k,u}^t is true.
        // We especially need t=0 at s, and some t=end at e.
        // First, find t_end such that D_{k,e}^t_end is true.
        int end_u = G.node(sc.lines[k].e.x, sc.lines[k].e.y);
        int start_u = G.node(sc.lines[k].s.x, sc.lines[k].s.y);

        // We don't know max t; search over the D-entries we actually have for this k.
        // Collect (t -> u) map for this k:
        unordered_map<int,int> t_to_u;

        // Iterate over all entries in vm.Dvar and filter by k
        for (const auto& it : vm.Dvar) {
            const auto& key = it.first;
            int var = it.second;
            if (key.k != k) continue;
            if (var >= (int)val.size()) continue;
            if (val[var] == 1) {
                // True
                auto got = t_to_u.find(key.t);
                // If multiple, prefer start/e matches
                if (got == t_to_u.end()) t_to_u.emplace(key.t, key.u);
            }
        }

        if (t_to_u.empty()) {
            // No info; return empty row -> will just print 0 below for UNSAT check outside
            paths[k] = {};
            continue;
        }

        // Find smallest t (expect 0) and t_end where u == end_u
        int t_min = INT_MAX, t_max = INT_MIN, t_end = -1;
        for (auto& p : t_to_u) {
            t_min = min(t_min, p.first);
            t_max = max(t_max, p.first);
            if (p.second == end_u) t_end = max(t_end, p.first);
        }
        if (t_min == INT_MAX) { paths[k] = {}; continue; }
        if (t_end < 0) t_end = t_max; // fallback if end wasn't found

        // Build directions from t_min..t_end
        vector<char> dirs;
        for (int t = t_min; t < t_end; ++t) {
            auto itA = t_to_u.find(t);
            auto itB = t_to_u.find(t+1);
            if (itA == t_to_u.end() || itB == t_to_u.end()) break;
            int u = itA->second, v = itB->second;
            auto [ux,uy] = G.xy(u);
            auto [vx,vy] = G.xy(v);
            if (vy == uy) {
                if (vx == ux+1) dirs.push_back('R');
                else if (vx == ux-1) dirs.push_back('L');
                else {
                    break;
                }
            } else if (vx == ux) {
                if (vy == uy+1) dirs.push_back('D');
                else if (vy == uy-1) dirs.push_back('U');
                else {
                    break;
                }
            } else {
                break;
            }
        }

        paths[k] = move(dirs);
    }

    return paths;
}

// ---------- Printing in the requested format ----------
static void printPaths(const vector<vector<char>>& paths, bool sat) {
    if (!sat) {
        cout << "0\n";
        return;
    }
    for (const auto& row : paths) {
        if (row.empty()) { cout << "0\n"; continue; }
        for (size_t i = 0; i < row.size(); ++i) {
            cout << row[i] << ' ';
        }
        cout << "0\n";
    }
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Args:
    //   argv[1] = model file (MiniSAT output), e.g. test.satoutput
    //   argv[2] = city file (default "test.city")
    //   argv[3] = varmap file (default "test.varmap")
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <test.satoutput> [test.city] [test.varmap]\n";
        return 1;
    }
    string modelPath = argv[1];
    string cityPath  = (argc >= 3 ? argv[2] : "test.city");
    string mapPath   = (argc >= 4 ? argv[3] : "test.varmap");

    // model
    vector<char> val; int maxVar = 0;
    bool sat = readModel(modelPath, val, maxVar);

    if (!sat) {
        cout << "0\n";
        return 0;
    }

    // Read scenario
    ifstream fin(cityPath);
    if (!fin) { cerr << "Cannot open " << cityPath << "\n"; return 1; }
    int scenarioType = 0;
    if (!(fin >> scenarioType)) { cerr << "Bad city file header\n"; return 1; }
    
    Scenario1 sc1;
    if (scenarioType == 1) {
        sc1 = parseScenario1(fin);
    } else if (scenarioType == 2) {
        Scenario2 sc2 = parseScenario2(fin);
        sc1.N = sc2.N; sc1.M = sc2.M; sc1.K = sc2.K; sc1.J = sc2.J;
        sc1.lines = sc2.lines;
    } else {
        cerr << "Unsupported scenario type (expected 1 or 2)\n";
        return 1;
    }

    Grid G(sc1.N, sc1.M);
    VarMap vm;
    vm.load(mapPath);
    auto paths = reconstructPaths(sc1, G, vm, val);

    printPaths(paths, true);
    return 0;
}