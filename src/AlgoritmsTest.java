package com.softwareTesting;
import java.util.*;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import GraphTesting.GraphAlgorithms;
import GraphTesting.GraphAlgorithms.Find_Centroid;
import GraphTesting.GraphAlgorithms.Graph;
import TreesTesting.TreesTraverse;
import TreesTesting.TreesTraverse.KruskalAlgorithm;
import TreesTesting.TreesTraverse.Trie;
import TreesTesting.TreesTraverseTest.BinaryTree;
import TreesTesting.TreesTraverseTest.Node;




class AlgoritmsTest {


	TreeAlgoritms TreesTester = new TreeAlgoritms();
	@Test	
    public void TestDFS1()
    {
        LinkedList<Integer> adj[];
        boolean visited[];  
        // make TreesTester case 
        int V = 5;
        adj = new LinkedList[V+1];  
        visited = new boolean[V+1];  
        
        for (int i = 0; i <= V; i++)  
            adj[i] = new LinkedList<Integer>(); 

        List<Integer> a = new ArrayList<Integer>();
        List<Integer> b = new ArrayList<Integer>();
        // edges b/w a[i] - b[i] for all i's
        int numEdges = 4;
        a.add(5);
        b.add(4);
        a.add(4);
        b.add(3);
        a.add(2);
        b.add(3);
        a.add(1);
        b.add(5);

        for(int i=0;i<numEdges;i++){
            adj[a.get(i)].add(b.get(i));
            adj[b.get(i)].add(a.get(i));
        }
        ArrayList<Integer> order = new ArrayList<Integer>();
        TreesTester.DFS(1, visited, order, adj);
        ArrayList<Integer> correct_order = new ArrayList<Integer>();
        correct_order.add(1);
        correct_order.add(5);
        correct_order.add(4);
        correct_order.add(3);
        correct_order.add(2);
        assertEquals(correct_order , order);    

    }
	@Test
    public void TestDFS2() {
		LinkedList<Integer> adj[];
        boolean visited[];  
        // make TreesTester case 
        int V = 8;
        adj = new LinkedList[V+1];  
        visited = new boolean[V+1];  
        
        for (int i = 0; i <= V; i++)  
            adj[i] = new LinkedList<Integer>(); 

        List<Integer> a = new ArrayList<Integer>();
        List<Integer> b = new ArrayList<Integer>();
        // edges b/w a[i] - b[i] for all i's
        int numEdges = 7;
        	a.add(5); b.add(2);
        	a.add(5); b.add(7);
        	a.add(1); b.add(7);
        	a.add(5); b.add(4);
        	a.add(3); b.add(8);
        	a.add(8); b.add(5);
        	a.add(7); b.add(6);
	
        for(int i=0;i<numEdges;i++){
            adj[a.get(i)].add(b.get(i));
            adj[b.get(i)].add(a.get(i));
        }
        ArrayList<Integer> order = new ArrayList<Integer>();
        TreesTester.DFS(1, visited, order, adj);
        ArrayList<Integer> correct = new ArrayList<Integer>();
        correct.add(1);
        correct.add(7);
        correct.add(5);
        correct.add(2);
        correct.add(4);
        correct.add(8);
        correct.add(3);
        correct.add(6);
        assertEquals(correct, order);        
      
	}
	
	@Test
	public void TestBFS1()
    {
        LinkedList<Integer> adj[];
        boolean visited[];  
        // make TreesTester case 
        int V = 5;
        adj = new LinkedList[V+1];  
        visited = new boolean[V+1];  
        
        for (int i = 0; i <= V; i++)  
            adj[i] = new LinkedList<Integer>(); 

        List<Integer> a = new ArrayList<Integer>();
        List<Integer> b = new ArrayList<Integer>();
        // edges b/w a[i] - b[i] for all i's
        int numEdges = 4;
        a.add(5);
        b.add(4);
        a.add(4);
        b.add(3);
        a.add(2);
        b.add(3);
        a.add(1);
        b.add(5);

        for(int i=0;i<numEdges;i++){
            adj[a.get(i)].add(b.get(i));
            adj[b.get(i)].add(a.get(i));
        }
        ArrayList<Integer> order = new ArrayList<Integer>();
        order = TreesTester.BFS(1, adj, V);
        ArrayList<Integer> correct_order = new ArrayList<Integer>();
        correct_order.add(1);
        correct_order.add(5);
        correct_order.add(4);
        correct_order.add(3);
        correct_order.add(2);
        assertEquals(correct_order , order);    

    }
	
	@Test
    public void TestBFS2() {
		LinkedList<Integer> adj[];
        boolean visited[];  
        // make TreesTester case 
        int V = 8;
        adj = new LinkedList[V+1];  
        visited = new boolean[V+1];  
        
        for (int i = 0; i <= V; i++)  
            adj[i] = new LinkedList<Integer>(); 

        List<Integer> a = new ArrayList<Integer>();
        List<Integer> b = new ArrayList<Integer>();
        // edges b/w a[i] - b[i] for all i's
        int numEdges = 7;
        	a.add(5); b.add(2);
        	a.add(5); b.add(7);
        	a.add(1); b.add(7);
        	a.add(5); b.add(4);
        	a.add(3); b.add(8);
        	a.add(8); b.add(5);
        	a.add(7); b.add(6);
	
        for(int i=0;i<numEdges;i++){
            adj[a.get(i)].add(b.get(i));
            adj[b.get(i)].add(a.get(i));
        }
        ArrayList<Integer> order = new ArrayList<Integer>();
        order = TreesTester.BFS(1, adj, V);
        ArrayList<Integer> correct = new ArrayList<Integer>();
        correct.add(1);
        correct.add(7);
        correct.add(5);
        correct.add(6);
        correct.add(2);
        correct.add(4);
        correct.add(8);
        correct.add(3);
        assertEquals(correct, order);        
      
	}
	

	@Test
    public void getMinWeight1(){
        int V = 5;
        int INF = Integer.MAX_VALUE;
        int graph[][] = new int[][] { { INF, 13, INF, 21, INF },
			{ 13, INF, 2, 4, 10 },
			{ INF, 2, INF, INF, 43 },
			{ 21, 4, INF, INF, 19 },
			{ INF, 10 ,43,19, INF } };
       assertEquals(2, TreesTester.getMinWeight(graph, V));
    }

	
	
	@Test
    public void getMinWeight2(){
        int V = 6;
        int INF = Integer.MAX_VALUE;
        int graph[][] = new int[][] { { INF, 9, INF, 4, INF, 2},
            {9, INF, 5, INF, INF, INF},
            {INF, 5, INF, INF, 3, INF},
            {4, INF, INF, INF, INF, INF},
            {INF, INF, 3, INF, INF, INF},
            {2, INF, INF, INF, INF, INF}};
       assertEquals(2, TreesTester.getMinWeight(graph, V));
    }


	@Test
    public void getMaxWeight1(){
        int V = 5;
        int INF = -Integer.MAX_VALUE;
        int graph[][] = new int[][] { { INF, 13, INF, 21, INF },
			{ 13, INF, 2, 4, 10 },
			{ INF, 2, INF, INF, 43 },
			{ 21, 4, INF, INF, 19 },
			{ INF, 10 ,43,19, INF } };
       assertEquals(43, TreesTester.getMaxWeight(graph, V));
    }

	
	
	@Test
    public void getMaxWeight2(){
        int V = 6;
        int INF = -Integer.MAX_VALUE;
        int graph[][] = new int[][] { { INF, 9, INF, 4, INF, 2},
            {9, INF, 5, INF, INF, INF},
            {INF, 5, INF, INF, 3, INF},
            {4, INF, INF, INF, INF, INF},
            {INF, INF, 3, INF, INF, INF},
            {2, INF, INF, INF, INF, INF}};
       assertEquals(9, TreesTester.getMaxWeight(graph, V));
    }
	
	
    class Node {
        int data;
        Node left, right;
     
        public Node(int item)
        {
            data = item;
            left = right = null;
        }
    }
    
    
    class BinaryTree {
        Node root;
     
        // Method to calculate the diameter and return it to
        // main
        int diameter(Node root)
        {
            // base case if tree is empty
            if (root == null)
                return 0;
     
            // get the height of left and right sub-trees
            int lheight = height(root.left);
            int rheight = height(root.right);
     
            // get the diameter of left and right sub-trees
            int ldiameter = diameter(root.left);
            int rdiameter = diameter(root.right);
     
            /* Return max of following three
              1) Diameter of left subtree
              2) Diameter of right subtree
              3) Height of left subtree + height of right
              subtree + 1
             */
            return Math.max(lheight + rheight + 1,
                            Math.max(ldiameter, rdiameter));
        }
     
        // A wrapper over diameter(Node root)
        int diameter() { return diameter(root); }
     
        // The function Compute the "height" of a tree. Height
        // is the number of nodes along the longest path from
        // the root node down to the farthest leaf node.
        static int height(Node node)
        {
            // base case tree is empty
            if (node == null)
                return 0;
     
            // If tree is not empty then height = 1 + max of
            //  left height and right heights
            return (1
                    + Math.max(height(node.left),
                               height(node.right)));
        }
    }
	
	@Test 
    public void TestDiameter1(){
        BinaryTree tree = new BinaryTree();
        tree.root = new Node(1);
        tree.root.left = new Node(2);
        tree.root.right = new Node(3);
        tree.root.left.left = new Node(4);
        tree.root.left.right = new Node(5);
        assertEquals(4, tree.diameter());
    }
	@Test 
    public void TestDiameter2(){
        BinaryTree tree = new BinaryTree();
        tree.root = new Node(1);
        tree.root.left = new Node(2);
        tree.root.left.left = new Node(3);
        tree.root.left.left.left = new Node(4);
        tree.root.left.left.left.left = new Node(5);
        assertEquals(5, tree.diameter());
    }
	
	@Test
	public void TestlargestRectangleArea()
    {

		int a1[] = {6, 2, 5, 4, 5, 1, 6}; 
		int a2[] = {3, 5, 1, 7, 5, 9};
		assertEquals(12,TreesTester.largestRectangleArea(a1));
		assertEquals(15,TreesTester.largestRectangleArea(a2));
		

    }

	
	@Test
	public void TestfindMedianSortedArrays1()
    {
		int a1[] = {-5, 3, 6, 12, 15};
		int a2[] = {-12, -10, -6, -3, 4, 10};
		int median = 3;
		int val = TreesTester.findMedianSortedArrays(a1, a2);
		assertEquals(median,val);

    }
	
	@Test
	public void TestfindMedianSortedArrays2()
    {
		int a1[] = {2, 3, 5, 8};
		int a2[] = {10, 12, 14, 16, 18, 20};
		int median = 11;
		int val = TreesTester.findMedianSortedArrays(a1, a2);
		assertEquals(median,val);

    }
	

	
	@Test
	public void TestKrushkalsAlgorithm() {
		int V = 6,E = 5;
		  com.softwareTesting.TreeAlgoritms.KruskalAlgorithm graph = TreesTester.new KruskalAlgorithm(V, E);
		  
	  
	       
	       graph.edgeArray[0].source = 0;
	       graph.edgeArray[0].destination = 1;
	       graph.edgeArray[0].weight = 9;
	       
	       graph.edgeArray[1].source = 0;
	       graph.edgeArray[1].destination = 3;
	       graph.edgeArray[1].weight = 4;

	       graph.edgeArray[2].source = 0;
	       graph.edgeArray[2].destination = 5;
	       graph.edgeArray[2].weight = 2;
	       
	       graph.edgeArray[3].source = 2;
	       graph.edgeArray[3].destination = 4;
	       graph.edgeArray[3].weight = 3;
	       
	       graph.edgeArray[4].source = 1;
	       graph.edgeArray[4].destination = 2;
	       graph.edgeArray[4].weight = 5;
	       
	            
	       
	       assertEquals(23, graph.applyKruskal());
	}
	
	@Test
	public void TestKrushkalsAlgorithm2() {
		int V = 5,E = 7;
		  com.softwareTesting.TreeAlgoritms.KruskalAlgorithm graph = TreesTester.new KruskalAlgorithm(V, E);
		  

 
	       graph.edgeArray[0].source = 0;
	       graph.edgeArray[0].destination = 1;
	       graph.edgeArray[0].weight = 13;
	       
	       graph.edgeArray[1].source = 1;
	       graph.edgeArray[1].destination = 2;
	       graph.edgeArray[1].weight = 2;

	       graph.edgeArray[2].source = 1;
	       graph.edgeArray[2].destination = 3;
	       graph.edgeArray[2].weight = 4;
	       
	       graph.edgeArray[3].source = 0;
	       graph.edgeArray[3].destination = 3;
	       graph.edgeArray[3].weight = 21;
	       
	       graph.edgeArray[4].source = 2;
	       graph.edgeArray[4].destination = 4;
	       graph.edgeArray[4].weight = 43;
	       
	       graph.edgeArray[5].source = 1;
	       graph.edgeArray[5].destination = 4;
	       graph.edgeArray[5].weight = 10;	            
	       
	       assertEquals(29, graph.applyKruskal());
	}
	
	
	@Test 
	public void TestTrie() {
		com.softwareTesting.TreeAlgoritms.Trie head = TreesTester.new Trie();
        head.insert("techie");
        head.insert("techi");
        head.insert("tech");
        
        assertEquals(true,head.search("tech"));
        assertEquals(true,head.search("techi"));
        assertEquals(true,head.search("techie"));
        assertEquals(false,head.search("techiedelight"));
        
        head.insert("techiedelight");
        assertEquals(true,head.search("tech"));
        assertEquals(true,head.search("techi"));
        assertEquals(true,head.search("techie"));
        assertEquals(true,head.search("techiedelight"));
        	
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		GraphAlgorithms GraphTester = new GraphAlgorithms();
		
		@Test
	    public void TestFlloyd(){
			int INF = 10000000;
			int graph[][] = { { 0, 21, INF, 32 },
	                          { INF, 0, 23, INF },
	                          { INF, INF, 0, 11 },
	                          { INF, INF, INF, 0 } };
	        int finalMatrix[][] = {{0, 21, 44, 32}, 
	                            {INF, 0, 23, 34}, 
	                            {INF, INF, 0, 11}, 
	                            {INF, INF, INF, 0}};
	        int V = 4;
	        assertArrayEquals(finalMatrix, GraphTester.floydWarshall(graph, V));
	    }
		
		@Test
	    public void TestPrims1(){
	        int V = 5;
	        int graph[][] = new int[][] { { 0, 13, 0, 21, 0 },
				{ 13, 0, 2, 4, 10 },
				{ 0, 2, 0, 0, 43 },
				{ 21, 4, 0, 0, 19 },
				{ 0, 10 ,43,19, 0 } };
	       assertEquals(29, GraphTester.primMST(graph, V));
	    }
	
		@Test
	    public void TestPrims2(){
	        int V = 6;
	        int graph[][] = new int[][] { { 0, 9, 0, 4, 0, 2},
	            {9, 0, 5, 0, 0, 0},
	            {0, 5, 0, 0, 3, 0},
	            {4, 0, 0, 0, 0, 0},
	            {0, 0, 3, 0, 0, 0},
	            {2, 0, 0, 0, 0, 0}};
	       assertEquals(23, GraphTester.primMST(graph, V));
	    }
		
	    @Test
	    public void TestDijsktra1(){
	        int V = 9;
	        int graph[][] = new int[][] { { 0, 4, 0, 0, 0, 0, 0, 8, 0 },
	                        { 4, 0, 8, 0, 0, 0, 0, 11, 0 },
	                        { 0, 8, 0, 7, 0, 4, 0, 0, 2 },
	                        { 0, 0, 7, 0, 9, 14, 0, 0, 0 },
	                        { 0, 0, 0, 9, 0, 10, 0, 0, 0 },
	                        { 0, 0, 4, 14, 10, 0, 2, 0, 0 },
	                        { 0, 0, 0, 0, 0, 2, 0, 1, 6 },
	                        { 8, 11, 0, 0, 0, 0, 1, 0, 7 },
	                        { 0, 0, 2, 0, 0, 0, 6, 7, 0 } };
	        int predicted[] = GraphTester.dijkstra(graph, 0, V);
	        int correct[] = {0, 4, 12, 19, 21, 11, 9, 8, 14};
	        assertArrayEquals(predicted, correct);
	    }
	    
	    @Test
	    public void TestDijsktra2(){
	        int V = 7;
	        int graph[][] = new int[][] { { 0, 0, 1, 2, 0, 0, 0 }, { 0, 0, 2, 0, 0, 3, 0 }, { 1, 2, 0, 1, 3, 0, 0 },
	            { 2, 0, 1, 0, 0, 0, 1 }, { 0, 0, 3, 0, 0, 2, 0 }, { 0, 3, 0, 0, 2, 0, 1 }, { 0, 0, 0, 1, 0, 1, 0 } };
	        int predicted[] = GraphTester.dijkstra(graph, 0, V);
	        int correct[] = {0,3,1,2,4,4,3};
	        assertArrayEquals(predicted, correct);
	    }
	    
	    @Test
	    public void testConnectedComp(){
	        ArrayList<ArrayList<Integer> > adjListArray;
	        adjListArray = new ArrayList<ArrayList<Integer>>();
	        int V = 5;
	        for (int i = 0; i < V; i++) {
				adjListArray.add(i, new ArrayList<>());
			}
	        // GraphTester case - 1
	        adjListArray.get(1).add(0);
	        adjListArray.get(0).add(1);
	        adjListArray.get(1).add(2);
	        adjListArray.get(2).add(1);
	        adjListArray.get(3).add(4);
	        adjListArray.get(4).add(3);
	        assertEquals(2, GraphTester.connectedComponents(adjListArray, V));
	
	
	    }
	    @Test
	    public void testConnectedComp2(){
	        ArrayList<ArrayList<Integer> > adjListArray;
	        adjListArray = new ArrayList<ArrayList<Integer>>();
	        int V = 8;
	        for (int i = 0; i < V; i++) {
				adjListArray.add(i, new ArrayList<>());
			}
		      adjListArray.get(1).add(2);
		      adjListArray.get(2).add(1);
		      adjListArray.get(3).add(4);
		      adjListArray.get(4).add(3);
		      adjListArray.get(5).add(6);
		      adjListArray.get(6).add(5);
		      adjListArray.get(7).add(6);
		      adjListArray.get(6).add(7);
		      assertEquals(4, GraphTester.connectedComponents(adjListArray, V));
	    }
	    
	
	
		
	    @Test 
	    public void TestisBc() {
	    	int V = 5;
	    	LinkedList<Integer> adj[];
	    	adj = new LinkedList[V];
	        for (int i=0; i<V; ++i)
	   	    	adj[i] = new LinkedList<Integer>();
	        
	       adj[1].add(0);
	       adj[0].add(1);
	       
	       adj[0].add(2);
	       adj[2].add(0);
	       
	       adj[2].add(1);
	       adj[1].add(2);
	       
	       adj[0].add(3);
	       adj[3].add(0);
	       
	       adj[3].add(4);
	       adj[4].add(3);
	
	       assertEquals(false,GraphTester.isBC(V,adj));
	    }
	    
	    
	    @Test 
	    public void TestisBc2() {
	    	int V = 5;
	    	LinkedList<Integer> adj[];
	    	adj = new LinkedList[V];
	        for (int i=0; i<V; ++i)
	   	    	adj[i] = new LinkedList<Integer>();
	        
	        
	        
	       adj[0].add(1);
	       adj[1].add(0);
	       
	       adj[1].add(2);
	       adj[2].add(1);
	       
	       adj[2].add(0);
	       adj[0].add(2);
	       
	       adj[0].add(3);
	       adj[3].add(0);
	       
	       adj[3].add(4);
	       adj[4].add(3);
	       
	       adj[2].add(4);
	       adj[4].add(2);
	       
	       assertEquals(true,GraphTester.isBC(V,adj));
	    }
	    
	    @Test
	    public void testMaxFlow(){
	       int V = 6;
	       int graph[][] = new int[][] {
	            { 0, 16, 13, 0, 0, 0 }, { 0, 0, 10, 12, 0, 0 },
	            { 0, 4, 0, 0, 14, 0 },  { 0, 0, 9, 0, 0, 20 },
	            { 0, 0, 0, 7, 0, 4 },   { 0, 0, 0, 0, 0, 0 }
	        };
	        assertEquals(23, GraphTester.fordFulkerson(graph, 0, 5, V));
	        assertEquals(19, GraphTester.fordFulkerson(graph, 0, 3, V));
	    }
	    
		@Test
		public void TestTopologicalSort1() {
			 ArrayList<ArrayList<Integer> > adj;
			 int V = 6;
			 adj = new ArrayList<ArrayList<Integer> >(V);
			 for (int i = 0; i < V; ++i)
		            adj.add(new ArrayList<Integer>());
			 adj.get(5).add(2);
			 adj.get(5).add(0);
			 adj.get(4).add(0);
			 adj.get(4).add(1);
			 adj.get(2).add(3);
			 adj.get(3).add(1);
			 assertArrayEquals(new int[]{5,4,2,3,1,0},GraphTester.topologicalSort(V,adj));   
		}
	    
		@Test
		public void TestTopologicalSort2() {
			 ArrayList<ArrayList<Integer> > adj;
			 int V = 7;
			 adj = new ArrayList<ArrayList<Integer> >(V);
			 for (int i = 0; i < V; ++i)
		            adj.add(new ArrayList<Integer>());
			 adj.get(5).add(2);
			 adj.get(5).add(0);
			 adj.get(4).add(0);
			 adj.get(4).add(1);
			 adj.get(2).add(3);
			 adj.get(3).add(1);
			 adj.get(6).add(3);
			 adj.get(6).add(0);
			 adj.get(5).add(3);
			 assertArrayEquals(new int[]{6,5,4,2,3,1,0},GraphTester.topologicalSort(V,adj));   
		}
		
		@Test 
	    public void TestBridge1(){
	        int V = 5;
	        LinkedList<Integer> adj[];
	        adj = new LinkedList[V];
	        for (int i=0; i<V; ++i)
	            adj[i] = new LinkedList();
	        adj[1].add(0);adj[0].add(1);
	        adj[2].add(0);adj[0].add(2);
	        adj[1].add(2);adj[2].add(1);
	        adj[3].add(0);adj[0].add(3);
	        adj[3].add(4);adj[4].add(3);
	        ArrayList<ArrayList<Integer>> bridges = new ArrayList<ArrayList<Integer>>();
	        GraphTester.bridge(adj, V, bridges);
	        ArrayList<ArrayList<Integer>> correct = new ArrayList<ArrayList<Integer>>();
	        ArrayList<Integer> b1 = new ArrayList<Integer>();        
	        b1.add(3);b1.add(4);
	        ArrayList<Integer> b2 = new ArrayList<Integer>();        
	        b2.add(0);b2.add(3);
	        correct.add(b1);
	        correct.add(b2);
	        assertEquals(bridges, correct);
	    }
		
		@Test 
	    public void TestBridge2(){
	        int V = 7;
	        LinkedList<Integer> adj[];
	        adj = new LinkedList[V];
	        for (int i=0; i<V; ++i)
	            adj[i] = new LinkedList();
	        adj[1].add(0);adj[0].add(1);
	        adj[2].add(1);adj[1].add(2);
	        adj[2].add(0);adj[0].add(2);
	        adj[3].add(1);adj[1].add(3);
	        adj[1].add(4);adj[4].add(1);
	        adj[1].add(6);adj[6].add(1);
	        adj[3].add(5);adj[5].add(3);
	        adj[4].add(5);adj[5].add(4);
	        
	        
	        ArrayList<ArrayList<Integer>> bridges = new ArrayList<ArrayList<Integer>>();
	        GraphTester.bridge(adj, V, bridges);
	        ArrayList<ArrayList<Integer>> correct = new ArrayList<ArrayList<Integer>>();
	        ArrayList<Integer> b1 = new ArrayList<Integer>();        
	        b1.add(1);
	        b1.add(6);
	        ArrayList<Integer> b2 = new ArrayList<Integer>();        
	        correct.add(b1);
	        assertEquals(bridges, correct);
	    }
		
		@Test 
	    public void TestSCC(){
	        int V = 5;
	        LinkedList<Integer> adj[];
	        adj = new LinkedList[V];
	        for (int i=0; i<V; ++i)
	            adj[i] = new LinkedList();
	        adj[1].add(0);
	        adj[0].add(2);
	        adj[2].add(1);
	        adj[0].add(3);
	        adj[3].add(4);
	        com.softwareTesting.GraphAlgorithms.Graph g = GraphTester.new Graph(V);
	        assertEquals(3, g.getSCCsCount(V, adj));
	    }
		
		
		@Test
		public void TestFindCentroid() {
	        int V = 4;
	                    // Number of nodes in the Tree
	        ArrayList<Integer>[] graph = new ArrayList[V];
	
	        for(int i=0;i<V;++i)
	            graph[i]=new ArrayList<>();
	                    //Initialisation
	
	
	
	        
	        graph[0].add(3);graph[3].add(0);
	        graph[3].add(2);graph[2].add(3);
	        graph[3].add(1);graph[1].add(3);
	        
	        
	        
	        com.softwareTesting.GraphAlgorithms.Find_Centroid fc = GraphTester.new Find_Centroid(graph,V);
	        
	        int centroid = fc.findCentroid(new java.util.Random().nextInt(V),graph);
	                    // Arbitrary indeed... ;)
	
	//        System.out.println("Centroid: "+(centroid+1));
	        assertEquals(centroid,3);
	                    // '+1' for output in 1-based index
		}


}
