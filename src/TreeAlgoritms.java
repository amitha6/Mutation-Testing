package com.softwareTesting;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import TreeAlgoritms.TreesTraverse.Node;
import TreeAlgoritms.TreesTraverse.Trie;
import TreeAlgoritms.TreesTraverse.KruskalAlgorithm.Edge;
import TreeAlgoritms.TreesTraverse.KruskalAlgorithm.Subset;

public class TreeAlgoritms {
	void DFS(int vertex, boolean vis[], ArrayList<Integer> order, LinkedList<Integer> adj[])
	{  
        
		vis[vertex] = true; 
        order.add(vertex);
        Iterator<Integer> it = adj[vertex].listIterator();  
        
        while (it.hasNext()) 
        {  
          int n = it.next();  
        
          if (!vis[n])
          {
        	  DFS(n, vis, order, adj);          	  
          }
        }
          
    }  




    ArrayList<Integer> BFS(int s, LinkedList<Integer> adj[],int V)
    {
        boolean visited[] = new boolean[V+1];
 
        LinkedList<Integer> queue = new LinkedList<Integer>();
        ArrayList<Integer> order = new ArrayList<Integer>();

        visited[s]=true;
        queue.add(s);
 
        while (queue.size() != 0)
        {

            s = queue.poll();

            order.add(s);
 
            Iterator<Integer> i = adj[s].listIterator();
            while (i.hasNext())
            {
                int n = i.next();
          
                if (!visited[n])
                {
                    visited[n] = true;
                    queue.add(n);
                }
            }
        }
        return order;
    }
    
    int getMinWeight(int graph[][], int V)
    {
        int minEdge = Integer.MAX_VALUE;
        for(int i=0;i<V;i++)
        {
            for(int j=0;j<V;j++)
            {
                minEdge = Math.min(minEdge, graph[i][j]);
            }
        }
        return minEdge; 
    }

    int getMaxWeight(int graph[][], int V)
    {
        int maxEdge = Integer.MIN_VALUE;
        
        for(int i=0;i<V;i++)
        {
            for(int j=0;j<V;j++)
            {
                maxEdge = Math.max(maxEdge, graph[i][j]);
            }
        }
        return maxEdge; 
    }
    
    class Node
    {
        int data;
        Node left, right;
     
        public Node(int item)
        {
            data = item;
            left = right = null;
        }
    }

    class BinaryTree 
    {
        Node root;
     
        int diameter(Node root)
        {
            if (root == null) 
            {
            	return 0;
            	
            }
     
            // get the height of left and right sub-trees
            int lheight = height(root.left);
            int rheight = height(root.right);
     
            // get the diameter of left and right sub-trees
            int ldiameter = diameter(root.left);
            int rdiameter = diameter(root.right);
     
            return Math.max(lheight + rheight + 1,
                            Math.max(ldiameter, rdiameter));
        }
     
        int diameter() 
        { 
        	return diameter(root); 
        }
     
        static int height(Node node)
        {
            if (node == null)
                return 0;
     
            return (1
                    + Math.max(height(node.left),
                               height(node.right)));
        }
    }
    
    
    public int largestRectangleArea(int[] heights) {
        int n = heights.length;
        if(n == 0) 
        {
        	return 0; 
        }
        int maxArea = 0;
        int left[] = new int[n]; 
        int right[] = new int[n]; 
        
        left[0] = -1;
        right[n - 1] = n;
        
        for(int i = 1; i < n; i++)
        {
            int prev = i - 1; 
            while(prev >= 0 && heights[prev] >= heights[i])
            {
                prev = left[prev]; 
            }
            left[i] = prev;
        }
      
        for(int i = n - 2; i >= 0; i--)
        {
            int prev = i + 1; 
            while(prev < n && heights[prev] >= heights[i])
            {
                prev = right[prev]; 
            }
            right[i] = prev;
        }
        
        for(int i = 0; i < n; i++)
        {
            int width = right[i] - left[i] - 1;
            maxArea = Math.max(maxArea, heights[i] * width);
        }
        return maxArea;
        
    }
    
    
	public int findMedianSortedArrays(int[] nums1, int[] nums2) {
	  int len1 = 0, len2 = 0;
	  
	  if(nums1 == null && nums2 == null)
	  {
	      return 0;
	  } 
	  else if(nums1 == null)
	  {
	      len2 = nums2.length;
	  } 
	  else if(nums2 == null) 
	  {
	      len1 = nums1.length;
	  } 
	  else
	  {
	      len1 = nums1.length;
	      len2 = nums2.length;
	  }
	  
	  if((len1 + len2) % 2 == 0) 
	  {
	      return (int)((findKth(nums1, 0, len1, nums2, 0, len2, (len1 + len2)/2) + findKth(nums1, 0, len1, nums2, 0, len2, (len1 + len2)/2 + 1))/2);
	  }
	  else
	  {
	      return (int)(findKth(nums1, 0, len1, nums2, 0, len2, (len1 + len2)/2 + 1));
	  }
	}

	public int findKth(int[] n1, int start1, int len1, int[] n2, int start2, int len2, int k) {
	  if(len1 > len2) 
	  {
		  return findKth(n2, start2, len2, n1, start1, len1, k);
	  }
	  
	  if(len1 == 0)
	  {
		  return n2[start2 + k - 1];
	  }
	  if(k == 1)
	  {
		  return Math.min(n1[start1], n2[start2]);
	  }
	  
	  int p1 = Math.min(len1, k / 2), p2 = k - p1;
	  int num1 = n1[start1 + p1 - 1], num2 = n2[start2 + p2 - 1];
	  
	  if(num1 == num2)
	  {
	      return num1;
	  } 
	  else if(num1 < num2) 
	  {
	      return findKth(n1, start1 + p1, len1 - p1, n2, start2, len2, k - p1);
	  } 
	  else 
	  {
	      return findKth(n1, start1, len1, n2, start2 + p2, len2 - p2, k - p2);
	  }
	}

    
    //Kruskals
	
	public class KruskalAlgorithm {  
	    class Edge implements Comparable<Edge>
	    {  
	        int source, destination, weight;  
	  
	        public int compareTo(Edge edgeToCompare) 
	        {  
	            return this.weight - edgeToCompare.weight;  
	        }  
	    };  
	  
	    class Subset 
	    {  
	        int parent, value;  
	    };  
	      
	    int vertices, edges;  
	    Edge edgeArray[];  
	  
	    KruskalAlgorithm(int vertices, int edges) 
	    {  
	        this.vertices = vertices;  
	        this.edges = edges;  
	        edgeArray = new Edge[this.edges];  
	        for (int i = 0; i < edges; ++i)  
	        {
	            edgeArray[i] = new Edge();  
	        }
	    }  
	      
	    public int applyKruskal() {  
	          
	        Edge finalResult[] = new Edge[vertices];  
	        int newEdge = 0;  
	        int index = 0;  
	        for (index = 0; index < vertices; ++index) 
	        {
	            finalResult[index] = new Edge(); 
	        }
	  
	        Arrays.sort(edgeArray);  
	          
	        Subset subsetArray[] = new Subset[vertices];  
	          
	        for (index = 0; index < vertices; ++index)
	        {
	            subsetArray[index] = new Subset();
	        }
	  
	        for (int vertex = 0; vertex < vertices; ++vertex) 
	        {  
	            subsetArray[vertex].parent = vertex;  
	            subsetArray[vertex].value = 0;  
	        }  
	        index = 0;  
	          
	        while (newEdge < vertices - 1) 
	        {  
	            Edge nextEdge = new Edge();  
	            nextEdge = edgeArray[index++];  
	              
	            int nextSource = findSetOfElement(subsetArray, nextEdge.source);  
	            int nextDestination = findSetOfElement(subsetArray, nextEdge.destination);  
	              
	            if (nextSource != nextDestination) 
	            {  
	                finalResult[newEdge++] = nextEdge;  
	                performUnion(subsetArray, nextSource, nextDestination);  
	            }  
	        }  
	        int answer = 0;
	        
	        for (index = 0; index < newEdge; ++index) 
	        {
	        	answer += finalResult[index].weight;
	        }
	        
	        return answer;
	        
	    }  
	      
	    int findSetOfElement(Subset subsetArray[], int i)
	    {  
	        if (subsetArray[i].parent != i) 
	        {  
	            subsetArray[i].parent = findSetOfElement(subsetArray, subsetArray[i].parent);  
	        }
	        return subsetArray[i].parent;  
	    }  
	      
	    void performUnion(Subset subsetArray[], int sourceRoot, int destinationRoot) {  
	          
	        int nextSourceRoot = findSetOfElement(subsetArray, sourceRoot);  
	        int nextDestinationRoot = findSetOfElement(subsetArray, destinationRoot);  
	  
	        if (subsetArray[nextSourceRoot].value < subsetArray[nextDestinationRoot].value)
	        {
	            subsetArray[nextSourceRoot].parent = nextDestinationRoot; 
	        }
	        else if (subsetArray[nextSourceRoot].value > subsetArray[nextDestinationRoot].value)
	        {
	            subsetArray[nextDestinationRoot].parent = nextSourceRoot; 
	        }
	        else 
	        {  
	            subsetArray[nextDestinationRoot].parent = nextSourceRoot;  
	            subsetArray[nextSourceRoot].value++;  
	        }  
	    }  
	}
	
	
	public class Trie
	{
	    private static final int CHAR_SIZE = 26;
	 
	    private boolean isLeaf;
	    private List<Trie> children = null;
	 
	    Trie()
	    {
	        isLeaf = false;
	        children = new ArrayList<>(Collections.nCopies(CHAR_SIZE, null));
	    }
	 
	    public void insert(String key)
	    {
	 
	        Trie curr = this;
	 
	        for (char c: key.toCharArray())
	        {
	            if (curr.children.get(c - 'a') == null)
	            {
	                curr.children.set(c - 'a', new Trie());
	            }
	 
	            curr = curr.children.get(c - 'a');
	        }
	 
	        curr.isLeaf = true;
	    }
	 
	    public boolean search(String key)
	    {
	 
	        Trie curr = this;
	 
	        for (char c: key.toCharArray())
	        {
	            curr = curr.children.get(c - 'a');
	 
	            if (curr == null)
	            {
	                return false;
	            }
	        }
	 
	        return curr.isLeaf;
	    }
	}
}
