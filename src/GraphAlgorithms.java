package com.softwareTesting;
import java.util.*;



public class GraphAlgorithms {
	int[][] floydWarshall(int dist[][], int V)
	{
        int i, j, k;
    
        for (k = 0; k < V; k++) 
        {
        	for (i = 0; i < V; i++)
        	{
            	for (j = 0; j < V; j++) 
            	{
            		if (dist[i][k] + dist[k][j] < dist[i][j])
            		{
            			dist[i][j] = dist[i][k] + dist[k][j];
            		}
            	}
        	}        
        }
    
        return dist;
    
	}
	

    int minKey(int key[], Boolean mstSet[], int V)
    {
        int min = Integer.MAX_VALUE, min_index = -1;
 
        for (int v = 0; v < V; v++)
        {
            if (mstSet[v] == false && key[v] < min) 
            {
                min = key[v];
                min_index = v;
            }
        }
 
        return min_index;
    }
 
    
    int totalMST(int key[], int V)
    {
        int sum = 0;
        for (int i = 0; i < V; i++)
        {
            sum += key[i];
        }
        return sum;
    }
 
    
    int primMST(int graph[][], int V)
    {
        int parent[] = new int[V];
 
        
        int key[] = new int[V];
 

        Boolean mstSet[] = new Boolean[V];
 
        for (int i = 0; i < V; i++) 
        {
            key[i] = Integer.MAX_VALUE;
            mstSet[i] = false;
        }
 

        key[0] = 0;

        parent[0] = -1; 
 
        
        for (int count = 0; count < V - 1; count++)
        {
            
            int u = minKey(key, mstSet, V);
 
            mstSet[u] = true;
 
            
            for (int v = 0; v < V; v++)
            {
 
            
                if (graph[u][v] != 0 && mstSet[v] == false && graph[u][v] < key[v]) 
                {
                    parent[v] = u;
                    key[v] = graph[u][v];
                }
            }
        }
 
        return totalMST(key, V);
    }

    int minDistance(int dist[], Boolean sptSet[], int V)
    {
        // Initialize min value
        int min = Integer.MAX_VALUE, min_index = -1;
 
        for (int v = 0; v < V; v++)
        {
            if (sptSet[v] == false && dist[v] <= min) 
            {
                min = dist[v];
                min_index = v;
            }
        }
        return min_index;
    }
    
    
    int[] dijkstra(int graph[][], int src, int V)
    {
        int dist[] = new int[V]; 
        Boolean sptSet[] = new Boolean[V];
 
        for (int i = 0; i < V; i++)
        {
            dist[i] = Integer.MAX_VALUE;
            sptSet[i] = false;
        }
 
        dist[src] = 0;
        for (int count = 0; count < V - 1; count++) 
        {
         
            int u = minDistance(dist, sptSet,V);
 
            
            sptSet[u] = true;
 
         
            for (int v = 0; v < V; v++)
                if (!sptSet[v] && graph[u][v] != 0
                    && dist[u] != Integer.MAX_VALUE
                    && dist[u] + graph[u][v] < dist[v])
                    dist[v] = dist[u] + graph[u][v];
        }
        return dist;
    }
    
    
    void DFSUtil(int v, boolean[] visited, ArrayList<ArrayList<Integer> > adjListArray)
	{
		// Mark the current node as visited and print it
		visited[v] = true;
		// Recur for all the vertices
		// adjacent to this vertex
		for (int x : adjListArray.get(v)) 
		{
			if (!visited[x])
			{
				DFSUtil(x, visited,adjListArray);
			}
		}
	}
	int connectedComponents(ArrayList<ArrayList<Integer> > adjListArray, int V)
	{
		// Mark all the vertices as not visited
		boolean[] visited = new boolean[V];
        int cnt =0 ;
        
		for (int v = 0; v < V; ++v)
		{
			if (!visited[v])
			{
				// print all reachable vertices
				// from v
				DFSUtil(v, visited,adjListArray);
				cnt += 1;
				// System.out.println();
			}
		}
		return cnt;
	}
	
	
	
	boolean isBCUtil(int u, boolean visited[], int disc[],int low[],int parent[], LinkedList<Integer> adj[],int time)
	{
		int children = 0;
		visited[u] = true;
		disc[u] = low[u] = ++time;
		Iterator<Integer> i = adj[u].iterator();
		
		
		while (i.hasNext())
		{
		   int v = i.next(); 
		   if (!visited[v])
		   {
		       children++;
		       parent[v] = u;
		       if (isBCUtil(v, visited, disc, low, parent,adj,time))
		       {
		    	   return true;
		       }
		       
		       low[u]  = Math.min(low[u], low[v]);
		       
		       if (parent[u] == -1 && children > 1)
		       {
		           return true;
		       }
		       if (parent[u] != -1 && low[v] >= disc[u])
		       {
		           return true;
		       }
		   }
		   else if (v != parent[u])
		   {
		       low[u]  = Math.min(low[u], disc[v]);
		   }
		}
		return false;
	}
	
	boolean isBC(int V, LinkedList<Integer> adj[])
	{
		boolean visited[] = new boolean[V];
		int disc[] = new int[V];
		int low[] = new int[V];
		int parent[] = new int[V];
		
		for (int i = 0; i < V; i++)
		{
		   parent[i] = -1;
		   visited[i] = false;
		}
		int time = 0;
		
		if (isBCUtil(0, visited, disc, low, parent, adj,time) == true)
		{
		   return false;
		}
		
		for (int i = 0; i < V; i++)
		{
		   if (visited[i] == false)
		   {
		       return false;
		   }
		}
		return true;
	}
	
	
	
	
	
	
	boolean bfs(int rGraph[][], int s, int t, int parent[], int V)
    {
        // Create a visited array and mark all vertices as
        // not visited
        boolean visited[] = new boolean[V];
        
        for (int i = 0; i < V; ++i)
        {
            visited[i] = false;
        }
        // Create a queue, enqueue source vertex and mark
        // source vertex as visited
        LinkedList<Integer> queue = new LinkedList<Integer>();
        
        queue.add(s);
        
        visited[s] = true;
        parent[s] = -1;
 
        // Standard BFS Loop
        while (queue.size() != 0) 
        {
            int u = queue.poll();
 
            for (int v = 0; v < V; v++) 
            {
                if (visited[v] == false && rGraph[u][v] > 0)
                {
                    if (v == t) 
                    {
                        parent[v] = u;
                        return true;
                    }
                    queue.add(v);
                    parent[v] = u;
                    visited[v] = true;
                }
            }
        }

        return false;
    }
 
    // Returns the maximum flow from s to t in the given
    // graph
    int fordFulkerson(int graph[][], int s, int t, int V)
    {
        int u, v;
 
        int rGraph[][] = new int[V][V];
 
        for (u = 0; u < V; u++)
        {
            for (v = 0; v < V; v++)
            {
                rGraph[u][v] = graph[u][v];
            }
        }
        // This array is filled by BFS and to store path
        int parent[] = new int[V];
 
        int max_flow = 0; // There is no flow initially
 
        // Augment the flow while there is path from source
        // to sink
        while (bfs(rGraph, s, t, parent, V)) {
 
            int path_flow = Integer.MAX_VALUE;
            
            for (v = t; v != s; v = parent[v]) 
            {
                u = parent[v];
                path_flow = Math.min(path_flow, rGraph[u][v]);
            }
 
            // update residual capacities of the edges and
            // reverse edges along the path
            for (v = t; v != s; v = parent[v]) 
            {
                u = parent[v];
                rGraph[u][v] -= path_flow;
                rGraph[v][u] += path_flow;
            }
 
            // Add path flow to overall flow
            max_flow += path_flow;
        }
 
        // Return the overall flow
        return max_flow;
    }
    
    void topologicalSortUtil(int v, boolean visited[],
            Stack<Integer> stack, ArrayList<ArrayList<Integer> > adj)
	{
		// Mark the current node as visited.
		visited[v] = true;
		Integer i;
		
		// Recur for all the vertices adjacent
		// to thisvertex
		Iterator<Integer> it = adj.get(v).iterator();
		while (it.hasNext()) 
		{
			i = it.next();
			if (!visited[i])
			{
				topologicalSortUtil(i, visited, stack, adj);
			}
		}
		
		// Push current vertex to stack
		// which stores result
		stack.push(v);
	}
	
	// The function to do Topological Sort.
	// It uses recursive topologicalSortUtil()
	int[] topologicalSort(int V, ArrayList<ArrayList<Integer> > adj)
	{
		Stack<Integer> stack = new Stack<Integer>();
		
		// Mark all the vertices as not visited
		boolean visited[] = new boolean[V];
		
		for (int i = 0; i < V; i++)
		{
			visited[i] = false;
		}
		
		// Call the recursive helper
		// function to store
		// Topological Sort starting
		// from all vertices one by one
		for (int i = 0; i < V; i++)
		{
			if (visited[i] == false)
			{
					topologicalSortUtil(i, visited, stack, adj);
			}
		}
		
		int []a = new int[stack.size()];

		int i = 0;
		
		while (stack.empty() == false) 
		{
			a[i++] = stack.pop();
		}
		
		return a;
	}
	
	void bridgeUtil(int u, boolean visited[], int disc[],int low[], int parent[], int time, ArrayList<ArrayList<Integer>> bridges,LinkedList<Integer> adj[])
	{
		
		// Mark the current node as visited
		visited[u] = true;
		
		// Initialize discovery time and low value
		disc[u] = low[u] = ++time;
		
		// Go through all vertices adjacent to this
		Iterator<Integer> i = adj[u].iterator();
		while (i.hasNext())
		{
		    int v = i.next();  // v is current adjacent of u
		
		    // If v is not visited yet, then make it a child
		    // of u in DFS tree and recur for it.
		    // If v is not visited yet, then recur for it
		    if (!visited[v])
		    {
		        parent[v] = u;
		        bridgeUtil(v, visited, disc, low, parent, time, bridges,adj);
		
		        // Check if the subtree rooted with v has a
		        // connection to one of the ancestors of u
		        low[u]  = Math.min(low[u], low[v]);
		
		        // If the lowest vertex reachable from subtree
		        // under v is below u in DFS tree, then u-v is
		        // a bridge
		        if (low[v] > disc[u])
		        {
		            ArrayList<Integer> inner = new ArrayList<Integer>();        
		            inner.add(u);
		            inner.add(v);
		            bridges.add(inner);
		        }
		    }
		
		    // Update low value of u for parent function calls.
		    else if (v != parent[u])
		    {
		        low[u]  = Math.min(low[u], disc[v]);
		    }
		}
	}


	void bridge(LinkedList<Integer> adj[], int V, ArrayList<ArrayList<Integer>> bridges)
	{
		boolean visited[] = new boolean[V];
		int disc[] = new int[V];
		int low[] = new int[V];
		int parent[] = new int[V];
		
		
		for (int i = 0; i < V; i++)
		{
		    parent[i] = -1;
		    visited[i] = false;
		}
		
		for (int i = 0; i < V; i++)
		{
		    if (visited[i] == false)
		    {
		        bridgeUtil(i, visited, disc, low, parent,0,bridges,adj);
		    }
		}
	}

	public class Graph
	{
	    private int V;   // No. of vertices
	    private LinkedList<Integer> adj[]; //Adjacency List
	 
	    //Constructor
	    Graph(int v)
	    {
	        V = v;
	        adj = new LinkedList[v];
	        for (int i=0; i<v; ++i)
	        {
	            adj[i] = new LinkedList();
	        }
	    }
	 
	    //Function to add an edge into the graph
//	    void addEdge(int v, int w)  { adj[v].add(w); }
	 
	    // A recursive function to print DFS starting from v
	    void DFSUtil(int v,boolean visited[],LinkedList<Integer> adj[])
	    {
	        // Mark the current node as visited and print it
	        visited[v] = true;
	 
	        int n;
	 
	        // Recur for all the vertices adjacent to this vertex
	        Iterator<Integer> i =adj[v].iterator();
	        while (i.hasNext())
	        {
	            n = i.next();
	            if (!visited[n])
	            {
	                DFSUtil(n,visited,adj);
	            }
	        }
	    }
	 
	    // Function that returns reverse (or transpose) of this graph
	    Graph getTranspose(LinkedList<Integer> adj[])
	    {
	        Graph g = new Graph(V);
	        for (int v = 0; v < V; v++)
	        {
	            // Recur for all the vertices adjacent to this vertex
	            Iterator<Integer> i =adj[v].listIterator();
	            while(i.hasNext())
	            {
	            	g.adj[i.next()].add(v);
	            }
	        }
	        return g;
	    }
	 
	    void fillOrder(int v, boolean visited[], Stack stack,LinkedList<Integer> adj[])
	    {
	        // Mark the current node as visited and print it
	        visited[v] = true;
	 
	        // Recur for all the vertices adjacent to this vertex
	        Iterator<Integer> i = adj[v].iterator();
	        while (i.hasNext())
	        {
	            int n = i.next();
	            if(!visited[n])
	            {
	                fillOrder(n, visited, stack,adj);
	            }
	        }
	 
	        stack.push(v);
	    }
	 
	    int getSCCsCount(int V,LinkedList<Integer> adj[])
	    {
	        Stack stack = new Stack();
	 
	        boolean visited[] = new boolean[V];
	        
	        for(int i = 0; i < V; i++)
	        {
	        	visited[i] = false;
	        }
	 
	        for (int i = 0; i < V; i++)
	        {
	            if (visited[i] == false) 
	            {
	                fillOrder(i, visited, stack,adj);
	            }
	        }
	        Graph gr = getTranspose(adj);
	 
	        for (int i = 0; i < V; i++)
	        {
	            visited[i] = false;
	        }
	 
	        int answer = 0;
	        
	        while (stack.empty() == false)
	        {
	            int v = (int)stack.pop();
	 
	            if (visited[v] == false)
	            {
	                gr.DFSUtil(v, visited,gr.adj);
	                answer++;
	            }
	        }
	        return answer;
	    }

	
	}
	
	
	public class Find_Centroid
	{
	    static final int MAXN=100_005;
	    ArrayList<Integer>[] graph;
	    static int[] depth,parent;  // Step 2
	    static int N;
	    
	    public Find_Centroid(ArrayList<Integer>[] graph,int N) 
	    {
	    	Find_Centroid.N = N;
	    	this.graph =  graph;
	    }

	   
	    static int[] queue=new int[MAXN],leftOver;
	                    // Step 3

	    static int findCentroid(int r,ArrayList<Integer>[] graph)
	    {
	        leftOver=new int[N];
	        int i,target=N/2,ach=-1;

	        bfs(r,graph);     // Step 4
	        for(i=N-1;i>=0;--i)
	            if(queue[i]!=r)
	                leftOver[parent[queue[i]]] += leftOver[queue[i]] +1;
	                    // Step 5
	        for(i=0;i<N;++i)
	            leftOver[i] = N-1 -leftOver[i];
	                    // Step 6
	        for(i=0;i<N;++i)
	            if(leftOver[i]<=target && leftOver[i]>ach)
	                    // Closest to target(=N/2) but does not exceed it.
	            {
	                r=i;    ach=leftOver[i];
	            }
	                    // Step 7
	        return r;
	    }
	    static void bfs(int root,ArrayList<Integer>[] graph)   // Iterative
	    {
	        parent=new int[N];  depth=new int[N];
	        int st=0,end=0;
	        parent[root]=-1;    depth[root]=1;
	                // Parent of root is obviously undefined. Hence -1.
	                // Assuming depth of root = 1
	        queue[end++]=root;
	        while(st<end)
	        {
	            int node = queue[st++], h = depth[node]+1;
	            Iterator<Integer> itr=graph[node].iterator();
	            while(itr.hasNext())
	            {
	                int ch=itr.next();
	                if(depth[ch]>0)     // 'ch' is parent of 'node'
	                    continue;
	                depth[ch]=h;   parent[ch]=node;
	                queue[end++]=ch;    // Recording the Traversal sequence
	            }
	        }
	    }
	}
	
 
}

