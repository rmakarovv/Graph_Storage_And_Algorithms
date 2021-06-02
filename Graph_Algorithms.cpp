/**
 * Graph ADT, Adjacency-Matrix Graph Implementation with the ability to
 * build graphs, transpose graphs, detect and print cycle in graphs,
 * find shortest path between two vertices of the graph.
 * Dijkstra's algorithm implementation.
 * @author Roman Makarov, github.com/rmakarovv
 */

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <queue>
#include <algorithm>

#define INF 10000007 // This constant is used for unknown path in Dijkstra algorithm
using namespace std;


/**
 * The class that represents a vertex in the graph.
 */
template <typename T>
class Vertex {
public:
    T name;
    string marked = "WHITE";

    Vertex() {};

    Vertex(Vertex *pVertex) {
        name = pVertex->name;
    }

    Vertex(Vertex const &vertex) {
        name = vertex.name;
    }

    Vertex(T s){
        name = s;
    }

    friend bool operator==(const Vertex& left, const Vertex& right){
        return left.name == right.name;
    }
};

/*
 * Operators for comparing vertices.
 */
template <typename T>
bool operator < (const Vertex<T> &v, const Vertex<T> &v2) {
    return ( v2.name < v.name );
}
template <typename T>
bool operator == (Vertex<T> &rhs, Vertex<T> &lhs) {
    return (rhs.name == lhs.name);
}

/**
 * The class that represents an Edge in the graph.
 * Inside the edge pointers to its ends are stored.
 * Also weight and bandwidth of an edge are used.
 */
template <typename T>
class Edge {
public:
    int weight;
    Vertex<T>* from;
    Vertex<T>* to;
    int bandwidth = 0;

    Edge(){}

    Edge(int w, Vertex<T>* f, Vertex<T>* t) {
        weight = w;
        from = f;
        to = t;
    }

    Edge(int w, Vertex<T>* f, Vertex<T>* t, int b) {
        weight = w;
        from = f;
        to = t;
        bandwidth = b;
    }

    Edge(Edge const &edge) {
        weight = edge.weight;
        from = edge.from;
        to = edge.to;
        bandwidth = edge.bandwidth;
    }

    friend bool operator==(const Edge& left, const Edge& right) {
        return (left.from == right.from && left.to == right.to && left.weight == right.weight);
    }
};

/**
 * The abstract class for Graph implementation
 * Inside general graph there are list of vertices and list of Edges.
 * All functions are virtual, therefore a class is abstract.
 */
template<typename T>
class Graph {
public:
    list<Vertex<T>> ver;
    list<Edge<T>> allEdges;

    Graph(){};

    virtual Vertex<T>& addVertex(T name) = 0;

    virtual void removeVertex(T name){};

    virtual Edge<T>& addEdge(int w, T f, T t) = 0;

    virtual Edge<T>& addEdgeBandwidth(int w, T f, T t, int b) = 0;

    virtual void removeEdge(T from, T to){};

    virtual vector<T> edgesFrom(T name){vector<T> v; return v;};

    virtual vector<Edge<T>> edgesTo(T name){vector<Edge<T>> v; return v;};

    virtual Vertex<T> findVertex(T name){return Vertex<T>();};

    virtual Edge<T> findEdge(T from, T to){return Edge<T>();}

    virtual bool vertexExist(T name){return true;};

    virtual bool hasEdge(T from, T to){return true;};
};


/**
 * The class that implements Adjacency Matrix Graph, extending
 * methods from "Graph" class.
 */
template <typename T>
class AdjacencyMatrixGraph : public Graph<T> {
public:
    // The representation of the adjacency matrix graph.
    // The names of the vertices are mapped to the information about them.
    // A pair of vertices is mapped to an edge between them if such edge exist.
    map<T, Vertex<T>> vertices;
    map<pair<T, T>, Edge<T>> edges;

    // Default constructor.
    AdjacencyMatrixGraph(){};

    // The function that adds vertex
    // with the check of already existing one.
    Vertex<T>& addVertex(T name) override {
        if (vertices.find(name) == vertices.end()) {
            Vertex<T> vertex(name);
            vertices[name] = vertex;
        }
        return vertices[name];
    }

    // The function that allows removing a vertex.
    // It also removes all edges that are connected to this vertex.
    void removeVertex(T name) override {
        vector<pair<T, T>> toDelete;
        for (const auto &el : edges) {
            if (el.first.first == name || el.first.second == name) {
                toDelete.push_back({el.first.first, el.first.second});
            }
        }
        vertices.erase(name);
        for (pair<T, T> &toDel : toDelete) {
            edges.erase({toDel.first, toDel.second});
        }
    }

    // The functions that allows adding edges.
    Edge<T>& addEdge(int w, T f, T t) override {
        Edge<T> edge(w, &vertices[f], &vertices[t]);
        edges[{f, t}] = edge;
        return edges[{f, t}];
    }
    // The functions that allows adding edges with specified bandwidth.
    Edge<T>& addEdgeBandwidth(int w, T f, T t, int b) override {
        Edge<T> edge(w, &vertices[f], &vertices[t], b);
        edges[{f, t}] = edge;
        return edges[{f, t}];
    }

    // The functions that allows removing edges that have a bandwidth
    // less than some specified value.
    void removeEdgeBandwidth(int b) {
        vector< pair<int,int> > toRemove;
        for (const auto &el : edges) {
            if (el.second.bandwidth < b)
                toRemove.push_back({el.first.first, el.first.second});
        }
        for (pair<int, int> &el : toRemove) {
            edges.erase(el);
        }
    }

    // The functions that allows removing any edge
    // if the endpoints are known.
    void removeEdge(T from, T to) override {
        edges.erase({from, to});
    }

    // The functions that returns all names of
    // the vertices that are adjacent to some vertex.
    vector<T> edgesFrom(T name) override {
        vector<T> ans;
        for (auto &i : edges) {
            if (i.first.first == name) {
                ans.push_back(i.first.second);
            }
        }
        return ans;
    }

    // The functions that returns all names of
    // the vertices that have a path to some vertex.
    vector<Edge<T>> edgesTo(T name) override {
        vector<Edge<T>> ans;
        for (auto &i : edges) {
            if (i.first.second == name) {
                ans.push_back(i.second);
            }
        }
        return ans;
    }

    // The functions that allows to know if there is an
    // edge between two vertices
    bool hasEdge(T from, T to) override {
        return (edges.find({from, to}) != edges.end());
    }

    // The function that transposes a graph.
    void transpose() {
        vector<pair<pair<T, T>, Edge<T>>> edgesToReverse;
        for (const auto &el : edges) {
            edgesToReverse.push_back({{el.first.first, el.first.second}, el.second});
        }

        edges.clear(); // clear the matrix to update it.

        for (int i = 0; i < edgesToReverse.size(); i++) {
            edges[{edgesToReverse[i].first.second, edgesToReverse[i].first.first}] = edgesToReverse[i].second;
        }
    }

    /**
     * The Depth First Search implemetation that returns true
     * if DFS is possible -> if the graph doesnt have cycles.
     *
     * start - starting vertex.
     * processed - stack of vertices to keep the cycle if it is found.
     *
     * If the vertex is marked "WHITE", then we did not yet considered it.
     * If the vertex is marked "GREY", then it is being processed now.
     * If the vertex is marked "BLACK" then is has already been processed.
     */
    bool DFS(Vertex<T> &start, stack<Vertex<T>> &processed) {
        processed.push(start);
        // If the vertex was already processed, then we are in the cycle
        if (start.marked == "GREY") {
            return false;
        }
        vector<T> adjacent = edgesFrom(start.name);

        start.marked = "GREY";

        // This is the way to detect if there is a cycle
        // by checking all adjacent vertices.
        for (T &ver : adjacent) {
            if (vertices[ver].marked == "GREY") {
                processed.push(vertices[ver]);
                return false;
            }
        }

        // Recursive DFS with the check for a return value.
        for (T &ver : adjacent) {
            bool ans = DFS(vertices[ver], processed);
            if ( !ans ) {
                return false;
            }
        }

        // Vertex that was already processed is marked black and deleted from stack
        start.marked = "BLACK";
        processed.pop();

        // If the DFS was successfully completed - true is returned.
        return true;
    }

    /**
     * This function allows to know if the graph contains a cycle.
     * If there are no cycle - it prints "ACYCLIC".
     * If some cycle is found the function prints its weight and involved vertices.
     */
    void isAcyclic() {
        stack<Vertex<T>> processed;
        bool isDFSPossible = true; // Detect if there is DFS possible a cycle
        // DFS is performed on each vertex to detect a cycle if it is present.
        for (auto &el : vertices) {
            // At each iteration vertices are marked "WHITE", and stack getting clear.
            for (auto &ver : vertices) {
                ver.second.marked = "WHITE";
            }
            while (!processed.empty())
                processed.pop();

            // DFS for each vertex.
            isDFSPossible = DFS(el.second, processed);
            // If cycle is found then the loop breaks out to print it.
            if (!isDFSPossible)
                break;
        }

        if (isDFSPossible) {
            cout << "ACYCLIC" << endl;
        } else {
            // Constructing a path of the loop
            Vertex<T> looped = processed.top(); // The vertex in which loop was detected.
            vector<Vertex<T>> path;
            // Pushing looped vertex into the path.
            path.push_back(processed.top());
            processed.pop();

            // While looped vertex is not found, path is filled.
            while (!(processed.top() == looped)) {
                path.push_back(processed.top());
                processed.pop();
            }
            path.push_back(path.front());
            // Path is reversed because previously is was in the stack
            reverse(path.begin(), path.end());

            // Weight is calculated by simply adding length of the edges.
            long long weight = 0;
            int i = 0;
            while ((i+1) != path.size()) {
                weight += edges[{path[i].name, path[i+1].name}].weight;
                i++;
            }
            path.pop_back();

            // Printing weight of the cycle and vertices inside it.
            cout << weight << " ";
            for (const Vertex<T>& ver : path) {
                cout << ver.name << " ";
            }
            cout << endl;
        }
    }
};


/**
 * This function is used for next commands:
 * 1) Add vertex / Add edge
 * 2) Remove vertex / Remove edge
 * 3) Check if there is an edge between two vertices (HAS_EDGE)
 * 4) Transpose a graph
 * 5) Check if the graph is acyclic.
 */
void general(){
    // Variables for command, names and weight.
    string command, name1, name2;
    int weight;

    // Initializing graph.
    AdjacencyMatrixGraph<string> graph;

    // Reading until the end of the file and performing commands
    while (cin >> command) {
        if (command == "ADD_VERTEX"){
            cin >> name1;
            graph.addVertex(name1);
        } else if (command == "ADD_EDGE") {
            cin >> name1 >> name2 >> weight;
            graph.addEdge(weight, name1, name2);
        } else if (command == "HAS_EDGE") {
            cin >> name1 >> name2;
            if (graph.hasEdge(name1, name2))
                cout << "TRUE" << endl;
            else cout << "FALSE" << endl;
        } else if (command == "REMOVE_EDGE") {
            cin >> name1 >> name2;
            graph.removeEdge(name1, name2);
        } else if (command == "REMOVE_VERTEX") {
            cin >> name1;
            graph.removeVertex(name1);
        } else if (command == "IS_ACYCLIC") {
            graph.isAcyclic();
        } else if (command == "TRANSPOSE") {
            graph.transpose();
        }
    }
}

/*
 * This is a comparator for a priority queue in which
 * first element is pair of {current Vertex,
 * and previousVertex(from which current was accessed).
 * The last parameter by which pairs are compared is length of the edge.
 */
class Compare{
public:
    bool operator() (pair< pair<int, int>,int> p1, pair< pair<int, int>,int> p2){
        return (p1.second > p2.second);
    }
};


/**
 * This is the function to find shortest path between two vertices
 * with respect to some minimum bandwidth
 * or print IMPOSSIBLE if there is no way between them
 * This function uses Dijkstra's algorithm.
 */
void shortestPath() {
    // Initializing number of vertices, number of queueOfVertices and graph.
    int n, m;
    cin >> n >> m;
    int source, target, length, bandwidth;
    AdjacencyMatrixGraph<int> graph;

    // Filling the graph with vertices
    for (int i = 1; i < n+1; i++) {
        graph.addVertex(i);
    }

    // Filling the graph with queueOfVertices
    for (int i = 0; i < m; i++) {
        cin >> source >> target >> length >> bandwidth;
        graph.addEdgeBandwidth(length, source, target, bandwidth);
    }

    // Initializing start and finish vertex, and minimal bandwidth.
    int start, finish, minBandwidth;
    cin >> start >> finish >> minBandwidth;

    // Removing all queueOfVertices that has bandwidth less than required
    graph.removeEdgeBandwidth(minBandwidth);

    // Initializing priority queue for vertices
    priority_queue<pair<pair<int, int>, int>, vector<pair<pair<int, int>, int>>, Compare > queueOfVertices;

    // The algorithm is build for vertices that have 'int' name.
    int current = start, previous = current; // current and previous vertices are assigned
    map<int, vector<int>> shortestPaths; // shortest path between source and each vertex
    map<int, int> shortestLength; // Shortest length between source and each vertex
    vector<int> path; // For creating a shortest path.

    // Initializing variables for starting vertex
    path.push_back(start);
    shortestPaths[start] = path;
    shortestLength[start] = 0;

    // To know if the vertex was already added to set
    vector<bool> marked(n+3, false);
    marked[start] = true;

    // All adjacent vertices of starting are pushed to the queue
    vector<int> adjacent = graph.edgesFrom(start);
    for (const int &el : adjacent) {
        queueOfVertices.push({{el, start}, graph.edges[{start, el}].weight});
    }

    // The loop that ends when a finish vertex is reached
    // or queue of available vertices is empty
    while (!queueOfVertices.empty() && !marked[finish]) {

        // Initializing previous and current vertex, and popping the queue.
        previous = queueOfVertices.top().first.second;
        current = queueOfVertices.top().first.first;
        queueOfVertices.pop();

        // If the vertex is was added to the set, proceed.
        if (!marked[current]) {
            // Updating shortest path between vertex and source.
            path = shortestPaths[previous]; path.push_back(current);
            shortestPaths[current] = path;
            // Updating shortest length between vertex and source.
            shortestLength[current] = shortestLength[previous] + graph.edges[{previous, current}].weight;
            marked[current] = true;

            // Getting all adjacent vertices of the current one
            adjacent = graph.edgesFrom(current);

            // We add adjacent vertices to the queue if they are not already in the set.
            for (int &el : adjacent) {
                if (!marked[el]) {
                    queueOfVertices.push({{el, current}, shortestLength[current] + graph.edges[{current, el}].weight });
                }
            }
        }

    }

    // If the finish vertex was reached.
    if (marked[finish]) {
        // Length between finish and start is simply stored in the created map.
        int lenght = shortestLength[finish];
        int bandwidthAns = -1; // initial value, assuming bandwidth is always positive
        path = shortestPaths[finish]; // Path is stored in the created map.
        int i = 0;
        // Calculating total bandwidth
        while ((i+1) != path.size()) {
            Edge<int> edge = graph.edges[{path[i], path[i+1]}];
            if (bandwidthAns == -1) {
                bandwidthAns = edge.bandwidth;
            } else {
                bandwidthAns = min(bandwidthAns, edge.bandwidth);
            }
            i++;
        }

        // Printing number of vertices in the path, its length and bandwidth
        cout << path.size() << " " << lenght << " " << bandwidthAns << endl;
        // Printing vertices in the path
        for (int &el : path) {
            cout << el << " ";
        }
        cout << endl;
    } else {
        // If there is no way between start and finish then print "IMPOSSIBLE".
        cout << "IMPOSSIBLE" << endl;
    }
}
