use std::collections::HashMap;
use std::f64::consts::PI;
use std::ffi::c_int;
use std::ffi::c_ulong;
use std::ffi::c_void;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::mem::size_of;
use std::ops::Index;

#[derive(Clone, Copy, Debug)]
struct Coord {
    field: [f64; 2],
}

impl Coord {
    fn lon(self) -> f64 {
        return self.field[0];
    }
    fn lat(self) -> f64 {
        return self.field[1];
    }
    fn with(self, i: usize, x: f64) -> Coord {
        assert!(i < 2);
        if i == 0 {
            return Coord {
                field: [x, self.lat()],
            };
        } else {
            return Coord {
                field: [self.lon(), x],
            };
        }
    }
    fn distance_to(self, other: Coord) -> f64 {
        let sin_dlon_2 = ((PI / 180f64) * (self.lon() - other.lon()) / 2f64).sin();
        let sin_dlat_2 = ((PI / 180f64) * (self.lat() - other.lat()) / 2f64).sin();
        let cos_lat1 = ((PI / 180f64) * self.lat()).cos();
        let cos_lat2 = ((PI / 180f64) * other.lat()).cos();
        let muladd = (sin_dlat_2 * sin_dlat_2) + (cos_lat1 * cos_lat2) * (sin_dlon_2 * sin_dlon_2);
        return 2f64 * muladd.sqrt().asin() * 6371000f64;
    }
}

impl Index<usize> for Coord {
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64 {
        return &self.field[i];
    }
}

#[derive(Debug)]
struct Record {
    field: (Coord, f64, i64),
}

impl Record {
    fn lon(&self) -> f64 {
        return self.field.0.lon();
    }
    fn lat(&self) -> f64 {
        return self.field.0.lat();
    }
    fn vel(&self) -> f64 {
        return self.field.1;
    }
    fn jid(&self) -> i64 {
        return self.field.2;
    }
    fn distance_to(&self, other: Coord) -> f64 {
        return self.field.0.distance_to(other);
    }
    fn coord(&self) -> Coord {
        return self.field.0;
    }
}
impl Index<usize> for Record {
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64 {
        return self.field.0.index(i);
    }
}

struct KdTreeNode {
    pub k: usize,
    pub x: f64,
    pub lo: Option<usize>,
    pub hi: Option<usize>,
}

struct KdTree<'a> {
    min_coord: Coord,
    max_coord: Coord,
    data: &'a Vec<Record>,
    nodes: Vec<KdTreeNode>,
    index: Vec<usize>,
}

trait NthElement<T: Sized, V: Copy> {
    fn nth_element<F>(&mut self, nth: T, less: &F)
    where
        F: Fn(&V, &V) -> bool;
}

impl<T: Copy> NthElement<usize, T> for &mut [T] {
    fn nth_element<F>(&mut self, nth: usize, less: &F)
    where
        F: Fn(&T, &T) -> bool,
    {
        if self.len() <= 1 {
            return;
        }
        // a   p   b
        let p = 0;
        let mut a = self.len() / 2;
        let mut b = self.len() - 1;
        if less(&self[p], &self[a]) {
            self.swap(p, a);
        }
        if less(&self[b], &self[p]) {
            self.swap(p, b);
        }
        if less(&self[p], &self[a]) {
            self.swap(p, a);
        }
        a = 0;
        b = self.len();
        while a + 1 < b {
            if less(&self[a + 1], &self[a]) {
                self.swap(a, a + 1);
                a += 1;
            } else {
                self.swap(a + 1, b - 1);
                b -= 1;
            }
        }
        if nth == a {
            return;
        }
        if nth < a {
            (&mut self[..a]).nth_element(nth, less);
        } else {
            (&mut self[b..]).nth_element(nth - b, less);
        }
    }
}

impl<'a> KdTree<'a> {
    fn variance(&self, k: usize, index: &[usize]) -> f64 {
        let n = index.len() as f64;
        let mut mean = 0 as f64;
        for i in index {
            mean += self.data[*i][k] / n;
        }
        let mut var = 0 as f64;
        for i in index {
            var += (self.data[*i][k] - mean) * (self.data[*i][k] - mean) / (n - 1f64);
        }
        return var;
    }
    fn add(&mut self, begin: usize, end: usize) -> Option<usize> {
        let index = &self.index[begin..end];
        if index.is_empty() {
            return None;
        }
        let n = index.len();
        let var0 = self.variance(0, index);
        let var1 = self.variance(1, index);
        let k = if var0 >= var1 { 0 } else { 1 };
        let mut index = &mut self.index[begin..end];
        index.nth_element(n / 2, &|a, b| self.data[*a][k] < self.data[*b][k]);
        index.swap(0, n / 2);
        let lon = self.data[index[0]].lon();
        let lat = self.data[index[0]].lat();
        if self.min_coord.lon().is_nan() || lon < self.min_coord.lon() {
            self.min_coord = self.min_coord.with(0, lon);
        }
        if self.max_coord.lon().is_nan() || lon > self.max_coord.lon() {
            self.max_coord = self.max_coord.with(0, lon);
        }
        if self.min_coord.lat().is_nan() || lat < self.min_coord.lat() {
            self.min_coord = self.min_coord.with(1, lat);
        }
        if self.max_coord.lat().is_nan() || lat > self.max_coord.lat() {
            self.max_coord = self.max_coord.with(1, lat);
        }
        self.nodes.push(KdTreeNode {
            k: k,
            x: self.data[index[0]][k],
            lo: None,
            hi: None,
        });
        let node_index = self.nodes.len() - 1;
        let lo = self.add(begin + 1, begin + n / 2 + 1);
        let hi = self.add(begin + n / 2 + 1, end);
        let node = &mut self.nodes[node_index];
        node.lo = lo;
        node.hi = hi;
        return Some(node_index);
    }
    fn new(data: &'a Vec<Record>) -> KdTree<'a> {
        let mut index: Vec<usize> = Vec::new();
        for i in 0..data.len() {
            index.push(i);
        }
        let mut tree = KdTree {
            min_coord: Coord {
                field: [f64::NAN, f64::NAN],
            },
            max_coord: Coord {
                field: [f64::NAN, f64::NAN],
            },
            data: data,
            nodes: Vec::new(),
            index: index,
        };
        tree.add(0, data.len());
        return tree;
    }
    fn range_<F>(
        &self,
        coord: Coord,
        node_index: usize,
        min_coord: Coord,
        max_coord: Coord,
        distance: f64,
        fun: &mut F,
    ) where
        F: FnMut(usize),
    {
        if self.data[self.index[node_index]].distance_to(coord) <= distance {
            fun(self.index[node_index]);
        }
        let node = &self.nodes[node_index];
        if let Some(lo) = node.lo {
            if !(coord[node.k] >= node.x
                && coord.distance_to(coord.with(node.k, node.x)) >= distance)
            {
                self.range_(
                    coord,
                    lo,
                    min_coord,
                    max_coord.with(node.k, node.x),
                    distance,
                    fun,
                );
            }
        }
        if let Some(hi) = node.hi {
            if !(coord[node.k] <= node.x
                && coord.distance_to(coord.with(node.k, node.x)) >= distance)
            {
                self.range_(
                    coord,
                    hi,
                    min_coord.with(node.k, node.x),
                    max_coord,
                    distance,
                    fun,
                );
            }
        }
    }
    fn range<F>(&self, coord: Coord, distance: f64, fun: &mut F)
    where
        F: FnMut(usize),
    {
        self.range_(coord, 0, self.min_coord, self.max_coord, distance, fun);
    }
}

fn read_data(filename: &str) -> Vec<Record> {
    let mut data: Vec<Record> = Vec::new();
    let file = File::open(filename).unwrap();
    for line in BufReader::new(file).lines() {
        if let Ok(content) = line {
            if !content.is_empty() {
                let record: Vec<&str> = content.split(',').collect();
                assert!(record.len() == 4);
                data.push(Record {
                    field: (
                        Coord {
                            field: [
                                record[0].parse::<f64>().unwrap(),
                                record[1].parse::<f64>().unwrap(),
                            ],
                        },
                        record[2].parse::<f64>().unwrap(),
                        record[3].parse::<f64>().unwrap() as i64,
                    ),
                });
                // println!("{:?}", data[data.len() - 1]);
            }
        }
    }
    return data;
}

fn filt_by_velocity(data: Vec<Record>) -> Vec<Record> {
    let n = data.len() as f64;
    let mut mean = 0 as f64;
    for record in &data {
        mean += record.vel() / n;
    }
    let mut var = 0 as f64;
    for record in &data {
        var += (record.vel() - mean) * (record.vel() - mean) / (n - 1f64);
    }
    let std = var.sqrt();
    let upper = mean + 3f64 * std;
    let mut data2: Vec<Record> = Vec::new();
    for record in data {
        if record.vel() <= upper {
            data2.push(record);
        }
    }
    return data2;
}

extern "C" {
    fn qsort(
        base: *mut c_void,
        nel: c_ulong,
        width: c_ulong,
        compar: extern "C" fn(*const c_void, *const c_void) -> c_int,
    );
}
static mut GRAPH: Vec<Vec<usize>> = Vec::new();
extern "C" fn compare(_pa: *const c_void, _pb: *const c_void) -> c_int {
    let a = unsafe { GRAPH[*(_pa as *const usize)].len() };
    let b = unsafe { GRAPH[*(_pb as *const usize)].len() };
    if a == b {
        return 0;
    }
    if a < b {
        return 1;
    } else {
        return -1;
    }
}

fn dfs(graph: &Vec<Vec<usize>>, minpts: usize, group: &mut Vec<i64>, i: usize, group_id: i64) {
    if group[i] != 0 {
        return;
    }
    group[i] = group_id;
    if graph[i].len() >= minpts {
        for j in &graph[i] {
            dfs(graph, minpts, group, *j, group_id);
        }
    }
}

fn main() {
    let data = read_data("../data.csv");
    let data = filt_by_velocity(data);
    let points = KdTree::new(&data);
    let eps = 8 as f64;
    let minpts = 50 as usize;
    let mut graph: Vec<Vec<usize>> = Vec::new();
    let mut min_distance: HashMap<i64, (usize, f64)> = HashMap::new();
    let mut pts: Vec<usize> = Vec::new();
    let mut core: Vec<usize> = Vec::new();
    for i in 0..data.len() {
        min_distance.clear();
        let coord = data[i].coord();
        points.range(coord, eps, &mut |j| {
            if data[i].jid() == data[j].jid() {
                return;
            }
            let distance = data[j].distance_to(coord);
            if let Some((index, min)) = min_distance.get_mut(&data[j].jid()) {
                if distance < *min || (distance == *min && j < *index) {
                    *index = j;
                    *min = distance;
                }
            } else {
                min_distance.insert(data[j].jid(), (j, distance));
            }
        });
        min_distance.insert(data[i].jid(), (i, 0 as f64));
        let mut dst: Vec<usize> = Vec::new();
        for entry in &min_distance {
            dst.push(entry.1 .0);
        }
        pts.push(dst.len());
        if dst.len() >= minpts {
            core.push(i);
        }
        graph.push(dst);
        if i % 5000usize == 0usize {
            println!("{}/{}", i, data.len());
        }
    }
    let graph = graph;
    let mut index: Vec<usize> = Vec::new();
    for i in 0..core.len() {
        index.push(i);
    }
    unsafe {
        GRAPH = graph.clone();
        qsort(
            (&mut index[..]).as_mut_ptr() as *mut c_void,
            index.len() as c_ulong,
            size_of::<usize>() as c_ulong,
            compare,
        )
    }
    let mut group = Vec::<i64>::new();
    let mut group_id = 0 as i64;
    for _ in 0..data.len() {
        group.push(0);
    }
    for idx in index {
        let i = core[idx];
        if graph[i].len() >= minpts && group[i] == 0 {
            group_id += 1;
            dfs(&graph, minpts, &mut group, i, group_id);
        }
    }
    let mut csv = File::create("../ans_rust.csv").unwrap();
    for i in 0..data.len() {
        writeln!(
            &mut csv,
            "{:.6},{:.6},{}",
            data[i].lon(),
            data[i].lat(),
            group[i]
        )
        .unwrap();
    }
    drop(csv);
}
