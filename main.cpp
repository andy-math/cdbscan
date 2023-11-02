#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numbers>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#ifdef assert
#undef assert
#endif
#define assert(condition, message)                                                                                     \
    {                                                                                                                  \
        if (!(condition)) {                                                                                            \
            myassert(__FILE__, __LINE__, condition, #condition, message);                                              \
        }                                                                                                              \
    }                                                                                                                  \
    (void)0
using namespace std;
void myassert(string file, int line, bool condition, string expr, string message) {
    cerr << "assertion fail at " << file << ':' << line << ": assert(" << expr << "); message: " << message << endl;
    exit(-1);
}

class Coord {
  public:
    double lon, lat;
    Coord(double lon, double lat) : lon(lon), lat(lat) {}
    double distanceTo(const Coord &other) const {
        double sin_dLat_2 = sin((numbers::pi / 180) * (lat - other.lat) / 2);
        double sin_dLon_2 = sin((numbers::pi / 180) * (lon - other.lon) / 2);
        double cos_lat1 = cos((numbers::pi / 180) * lat);
        double cos_lat2 = cos((numbers::pi / 180) * other.lat);
        double muladd = (sin_dLat_2 * sin_dLat_2) + (cos_lat1 * cos_lat2) * (sin_dLon_2 * sin_dLon_2);
        return 2 * asin(sqrt(muladd)) * 6371000;
    }
    double operator[](int i) const { return i == 0 ? lon : lat; }
    Coord with(int i, double x) { return i == 0 ? Coord{x, lat} : Coord{lon, x}; };
};

class Record : public Coord {
  public:
    const double vel, jid;
    Record(double lon, double lat, double vel, double jid) : Coord(lon, lat), vel(vel), jid(jid) {}
};

template <> class vector<Record> {
    vector<double> lon;
    vector<double> lat;
    vector<double> vel;
    vector<double> jid;
    vector<Record>(const vector<Record> &) = delete;
    vector<Record> &operator=(const vector<Record> &) = delete;
    vector<Record> &operator=(vector<Record> &&) = delete;

  public:
    vector<Record>() = default;
    vector<Record>(vector<Record> &&) = default;
    size_t size() const { return lon.size(); }
    void emplace_back(double lon, double lat, double vel, double jid) {
        this->lon.emplace_back(lon);
        this->lat.emplace_back(lat);
        this->vel.emplace_back(vel);
        this->jid.emplace_back(jid);
    }
    Record operator[](int i) const { return Record(lon[i], lat[i], vel[i], jid[i]); }
};

class KdTreeNode {
  public:
    const int k;
    const double x;
    int lo;
    int hi;
    KdTreeNode(int k, double x, int lo, int hi) : k(k), x(x), lo(lo), hi(hi) {}
};

template <> class vector<KdTreeNode> {
    vector<int> k;
    vector<double> x;
    vector<int> lo;
    vector<int> hi;
    vector<KdTreeNode>(const vector<KdTreeNode> &) = delete;
    vector<KdTreeNode> &operator=(const vector<KdTreeNode> &) = delete;
    vector<KdTreeNode> &operator=(vector<KdTreeNode> &&) = delete;

  public:
    vector<KdTreeNode>(vector<KdTreeNode> &&) = default;
    vector<KdTreeNode>() = default;
    size_t size() { return k.size(); }
    KdTreeNode operator[](int i) const { return KdTreeNode(k[i], x[i], lo[i], hi[i]); }
    void emplace_back(int k, double x) {
        this->k.emplace_back(k);
        this->x.emplace_back(x);
        this->lo.emplace_back(-1);
        this->hi.emplace_back(-1);
    }
    void set(int i, KdTreeNode node) {
        lo[i] = node.lo;
        hi[i] = node.hi;
    }
};

class KdTree {
    vector<KdTreeNode> nodes;
    vector<int> index;
    double minLon, maxLon, minLat, maxLat;
    double variance(const vector<Record> &data, int k, vector<int>::iterator begin, vector<int>::iterator end) {
        const double n = end - begin;
        double mean = 0;
        for (auto it = begin; it < end; it++) {
            mean += data[*it][k] / n;
        }
        double var = 0;
        for (auto it = begin; it < end; it++) {
            var += (data[*it][k] - mean) * (data[*it][k] - mean) / n;
        }
        return var;
    }

    int add(const vector<Record> &data, vector<int>::iterator begin, vector<int>::iterator end) {
        if (begin == end) {
            return -1;
        }
        const double var0 = variance(data, 0, begin, end);
        const double var1 = variance(data, 1, begin, end);
        const int k = var0 >= var1 ? 0 : 1;
        nth_element(begin, begin + (end - begin) / 2, end,
                    [k, &data](auto i, auto j) { return data[i][k] < data[j][k]; });
        auto tmp = *begin;
        *begin = *(begin + (end - begin) / 2);
        *(begin + (end - begin) / 2) = tmp;
        minLon = min(minLon, data[*begin].lon);
        maxLon = max(maxLon, data[*begin].lon);
        minLat = min(minLat, data[*begin].lat);
        maxLat = max(maxLat, data[*begin].lat);
        nodes.emplace_back(k, data[*begin][k]);
        int index = nodes.size() - 1;
        KdTreeNode node = nodes[index];
        node.lo = add(data, begin + 1, begin + (end - begin) / 2 + 1);
        node.hi = add(data, begin + (end - begin) / 2 + 1, end);
        nodes.set(index, node);
        return index;
    }
    template <typename T>
    void range(const vector<Record> &data, Coord coord, int nodeIndex, Coord minCoord, Coord maxCoord, double distance,
               const T &fun) const {
        assert(nodeIndex >= 0, "nodeIndex must greater than zero");
        if (data[index[nodeIndex]].distanceTo(coord) <= distance) {
            fun(index[nodeIndex]);
        }
        KdTreeNode node = nodes[nodeIndex];
        if (node.lo != -1 && !(coord[node.k] >= node.x && coord.distanceTo(coord.with(node.k, node.x)) >= distance)) {
            range(data, coord, node.lo, minCoord, maxCoord.with(node.k, node.x), distance, fun);
        }
        if (node.hi != -1 && !(coord[node.k] <= node.x && coord.distanceTo(coord.with(node.k, node.x)) >= distance)) {
            range(data, coord, node.hi, minCoord.with(node.k, node.x), maxCoord, distance, fun);
        }
    }

  public:
    KdTree(const vector<Record> &data)
        : minLon(data[0].lon), maxLon(data[0].lon), minLat(data[0].lat), maxLat(data[0].lat) {
        index.reserve(data.size());
        for (int i = 0; i < data.size(); i++) {
            index.emplace_back(i);
        }
        add(data, index.begin(), index.end());
    }
    template <typename T> void range(const vector<Record> &data, Coord coord, double distance, const T &fun) const {
        range(data, coord, 0, Coord{minLon, minLat}, Coord{maxLon, maxLat}, distance, fun);
    }
};

vector<Record> readData(string filename) {
    ifstream file{filename};
    string content{istreambuf_iterator<char>(file), istreambuf_iterator<char>()};
    vector<Record> data;
    for (const auto it : views::split(content, string_view("\n"))) {
        string_view line{it.begin(), it.end()};
        if (!line.empty()) {
            int i = 0;
            array<double, 4> buf;
            for (const auto x : views::split(line, string_view(","))) {
                assert(i < 4, "line constains incorrect numbers");
                buf[i++] = stod(string({x.begin(), x.end()}));
            }
            assert(i == 4, "line constains incorrect numbers");
            data.emplace_back(buf[0], buf[1], buf[2], buf[3]);
        }
    }
    for (int i = 1; i < data.size(); i++) {
        assert(data[i - 1].jid <= data[i].jid, "jid must be sorted");
    }
    return data;
}

vector<Record> filtByVelocity(const vector<Record> &data) {
    size_t n = data.size();
    double mean = 0;
    for (int i = 0; i < n; i++) {
        mean += data[i].vel / n;
    }
    double var = 0;
    for (int i = 0; i < n; i++) {
        var += (data[i].vel - mean) * (data[i].vel - mean) / (n - 1);
    }
    double std = sqrt(var);
    double upper = mean + 3 * std;
    vector<Record> data2;
    for (int i = 0; i < n; i++) {
        if (data[i].vel <= upper) {
            data2.emplace_back(data[i].lon, data[i].lat, data[i].vel, data[i].jid);
        }
    }
    return data2;
}

unordered_map<int, vector<int>> graph;

void dfs(const vector<Record> &data, int minpts, vector<int> &group, int i, int groupId) {
    if (group[i] != 0) {
        return;
    }
    group[i] = groupId;
    if (graph.at(i).size() >= minpts) {
        for (int j : graph.at(i)) {
            dfs(data, minpts, group, j, groupId);
        }
    }
}

int main() {
    vector<Record> data{filtByVelocity(readData("data.csv"))};
    KdTree points{data};
    double eps = 8;
    int minpts = 50;
    vector<int> core;
    unordered_map<double, pair<int, double>> minDistance;
    for (int i = 0; i < data.size(); i++) {
        minDistance.clear();
        points.range(data, data[i], eps, [&minDistance, &data, i](int j) {
            if (data[i].jid == data[j].jid) {
                return;
            }
            double distance{data[i].distanceTo(data[j])};
            if (minDistance.contains(data[j].jid)) {
                auto [index, min] = minDistance[data[j].jid];
                if (distance < min || (distance == min && j < index)) {
                    minDistance[data[j].jid] = {j, distance};
                }
            } else {
                minDistance[data[j].jid] = {j, distance};
            }
        });
        minDistance.emplace(data[i].jid, pair<int, double>{i, 0});
        for (auto &&[key, value] : minDistance) {
            graph[i].emplace_back(value.first);
        }
        int count = graph[i].size();
        if (count >= minpts) {
            core.emplace_back(i);
        }
        if (i % 5000 == 0) {
            cout << i << '/' << data.size() << endl;
        }
    }
    vector<int> index;
    index.reserve(data.size());
    for (int i = 0; i < core.size(); i++) {
        index.emplace_back(i);
    }
    qsort(&index[0], index.size(), sizeof(int), [](const void *pa, const void *pb) {
        size_t a = graph[*static_cast<const int *>(pa)].size();
        size_t b = graph[*static_cast<const int *>(pb)].size();
        return a == b ? 0 : (a < b ? 1 : -1);
    });
    vector<int> group;
    group.reserve(data.size());
    for (int i = 0; i < data.size(); i++) {
        group.emplace_back(0);
    }
    int groupId{0};
    for (int idx : index) {
        int i = core[idx];
        if (graph[i].size() >= minpts && !group[i]) {
            groupId++;
            dfs(data, minpts, group, i, groupId);
        }
    }
    FILE *csv = fopen("ans_compare.csv", "w");
    for (int i = 0; i < data.size(); i++) {
        fprintf(csv, "%lf,%lf,%d\n", data[i].lon, data[i].lat, group[i]);
    }
    fclose(csv);
    return 0;
}