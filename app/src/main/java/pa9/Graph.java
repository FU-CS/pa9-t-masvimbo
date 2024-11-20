package pa9;
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
    public int[] minimumSpanningTree();

    /**
     * Find the minimum spanning tree of the graph, using Prim's algorithm.
     */
    public int[] minimumSpanningTreePrim();

}
