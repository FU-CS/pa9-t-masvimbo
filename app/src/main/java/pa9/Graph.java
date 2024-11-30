package pa9;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
// import java.util.LinkedList;
// import java.util.Queue;

/** 
 * Graph interface for a directed, connected graph with vertices numbered from 0 to n-1.
 * For weighted graphs, the default weight of each edge is assumed to be 1.
 */
interface Graph {

    /**
     *  Add a weighted edge between two vertices.
     *  @param v vertex 1 (0-indexed)
     *  @param w vertex 2 (0-indexed)
     *  @param weight the weight of the edge
     * 
     * When called for unweighted graphs, it should be equivalent to addEdge(v, w).
     * 
     */
    public void addWeightedEdge(int v, int w, int weight);

    /**
     * Determine if there is a cycle in the graph.
     * @return true if the graph contains a cycle, false otherwise
     */ 
    public boolean hasNegativeCycle();

    /**
     * Find the shortest path between v and all other nodes, using Bellman Ford. 
     * 
     * @param  v source vertex
     * @return an array of integers representing the minimum cost 
     *         of the shortest path from v to all other vertices
     */ 
    public int[] shortestPath(int v);

    /**
     * Find the minimum spanning tree of the graph, using Kruskal's algorithm.
     */
    public HashSet<Edge> minimumSpanningTree();

    /**
     * Find the minimum spanning tree of the graph, using Prim's algorithm.
     */
    public int[] minimumSpanningTreePrim();



    public static class Edge implements Comparable<Edge> {
        int source;
        int destination;
        int weight;

        public Edge(int source, int destination, int weight) {
            this.source = source;
            this.destination = destination;
            this.weight = weight;
        }

        @Override
        public int compareTo(Edge other) {
            return Integer.compare(this.weight, other.weight);
        }
    }

    class GraphAdjacencyList implements Graph {
        private int vertices;
        private List<List<Edge>> adjacencyList;

        public GraphAdjacencyList(int vertices) {
            this.vertices = vertices;
            adjacencyList = new ArrayList<>(vertices);
            for (int i = 0; i < vertices; i++) {
                adjacencyList.add(new ArrayList<>());
            }
        }

        @Override
        public void addWeightedEdge(int v, int w, int weight) {
            adjacencyList.get(v).add(new Edge(v, w, weight));
        }

    @Override
    public boolean hasNegativeCycle() {
    int[] distances = new int[vertices];
    Arrays.fill(distances, Integer.MAX_VALUE);
    distances[0] = 0; // Default to vertex 0 as the starting point

    // Relax all edges |V| - 1 times
    for (int i = 0; i < vertices - 1; i++) {
        for (int u = 0; u < vertices; u++) {
            for (Edge edge : adjacencyList.get(u)) {
                int newDist = distances[u] + edge.weight;
                if (distances[u] != Integer.MAX_VALUE && newDist < distances[edge.destination]) {
                    distances[edge.destination] = newDist;
                }
            }
        }
    }
    for (int u = 0; u < vertices; u++) {
        for (Edge edge : adjacencyList.get(u)) {
            if (distances[u] != Integer.MAX_VALUE && distances[u] + edge.weight < distances[edge.destination]) {
                return true;
            }
        }
    }
    return false;
}

    
    @Override
    public int[] shortestPath(int v) {
    if (this.vertices == 0) {
        return new int[] {0};
    }
    int[] distances = new int[vertices];
    Arrays.fill(distances, Integer.MAX_VALUE);
    distances[v] = 0;

    // Relax all edges |V| - 1 times
    for (int i = 0; i < vertices - 1; i++) {
        // i is looping over every inner list in the adjacency list
        for (int u = 0; u < vertices; u++) {
            // looping over evey edge
            for (Edge edge : adjacencyList.get(u)) {
                int newDist = distances[u] + edge.weight;
                if (distances[u] != Integer.MAX_VALUE && newDist < distances[edge.destination]) {
                    distances[edge.destination] = newDist;
                }
            }
        }
    }
    return distances;
    }

        
    @Override
    public HashSet<Edge> minimumSpanningTree() {
    // PriorityQueue to store edges, automatically sorted by weight
    if (this.vertices == 0) {
        HashSet<Edge> defaultSet = new HashSet<>();
        defaultSet.add(new Edge(0, 0, 0));
        return defaultSet;
    }
    PriorityQueue<Edge> pq = new PriorityQueue<>();
    for (int u = 0; u < vertices; u++) {
        for (Edge edge : adjacencyList.get(u)) {
            pq.add(edge);
        }
    }

    // Initialize an empty graph to store the MST
    GraphAdjacencyList mst = new GraphAdjacencyList(vertices);
    HashSet<Edge> mstEdges = new HashSet<>();

    // Kruskal's Algorithm
    while (!pq.isEmpty() && mstEdges.size() < vertices - 1) {
        Edge edge = pq.poll();
        
        // Add edge to MST graph and check for cycles
        mst.addWeightedEdge(edge.source, edge.destination, edge.weight);
        if (mst.hasNegativeCycle()) {
            // Remove the edge if it creates a cycle
            mst.adjacencyList.get(edge.source).removeIf(e -> e.destination == edge.destination && e.weight == edge.weight);
        } else {
            // Otherwise, add it to the MST edges list
            mstEdges.add(edge);
        }
    }

    // If the MST is incomplete, return an empty set
    if (mstEdges.size() != vertices - 1) {
        return new HashSet<>(); // No MST exists
    }

    return mstEdges;
}

    @Override
    public int[] minimumSpanningTreePrim() {
        if (this.vertices == 0) {
            return new int[] {0};
        }
        // Prim - Jarnik's Algorithm
        int[] parent = new int[vertices];
        int[] distances = new int[vertices];
        boolean[] inMST = new boolean[vertices];
    
        Arrays.fill(distances, Integer.MAX_VALUE);
        distances[0] = 0;
        parent[0] = -1;
    
        PriorityQueue<Edge> pq = new PriorityQueue<>();
        pq.add(new Edge(-1, 0, 0));
    
        while (!pq.isEmpty()) {
            Edge current = pq.poll();
            int u = current.destination;
    
            if (inMST[u]) continue;
    
            inMST[u] = true;
            if (current.source != -1) {
                parent[u] = current.source;
            }
    
            for (Edge edge : adjacencyList.get(u)) {
                int v = edge.destination;
                int weight = edge.weight;
    
                if (!inMST[v] && weight < distances[v]) {
                    distances[v] = weight;
                    pq.add(new Edge(u, v, weight));
                }
            }
        }
        return parent;
    }        
}
public class Main {
    public static void main(String[] args) {
        Graph.GraphAdjacencyList graph = new Graph.GraphAdjacencyList(5);

        graph.addWeightedEdge(0, 1, 2);
        graph.addWeightedEdge(0, 2, 4);
        graph.addWeightedEdge(1, 2, 1);
        graph.addWeightedEdge(1, 3, 7);
        graph.addWeightedEdge(2, 4, 3);

        int[] shortestPath = graph.shortestPath(0);
        for (int distance : shortestPath) {
            System.out.print(distance + " ");
        }
        System.out.println();

        boolean hasNegativeCycle = graph.hasNegativeCycle();
        System.out.println(hasNegativeCycle);

        HashSet<Graph.Edge> mstEdges = graph.minimumSpanningTree();
        for (Graph.Edge edge : mstEdges) {
            System.out.println(edge.source + " -> " + edge.destination + " : " + edge.weight);
        }

        int[] mstPrim = graph.minimumSpanningTreePrim();
        for (int i = 0; i < mstPrim.length; i++) {
            System.out.println(i + " <- " + mstPrim[i]);
        }
    }
}
}