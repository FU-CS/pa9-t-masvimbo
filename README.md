[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/ZNKI4Thm)


<center><h1>PROGRAMMING ASSIGNMENT 9
<br/><br/>
Shortest Path and Minimum Spanning Tree</h1>
<br/>

**Due: Wednesday, November 27, 2024, 2:30 PM**</center>

<br/>

In this assignment, you are to implement the `Graph` interface in the `Graph.java` file provided in `app/src/main/java/pa9` folder. The interface defines the methods that you need to implement. 

<u>Your implementation should be based on the <b>Adjacency List</b> representation of a graph</u>.

Make sure to read the java doc comments in the `Graph` interface to understand the requirements of each method.

In total, you are to implement four methods in the `Graph` interface. 

--- 

## Q1. Bellman Ford Algorithm

### 1a. Shortest Path 

In this assignment, you are to implement the [Bellman Ford algorithm](https://fahadsultan.com/csc223/datastructs/graphs_bellmanford.html) to find the shortest path from a source vertex to all other vertices in the graph.

### 1b. Negative Cycle Detection

In addition to finding the shortest path, you are to implement a method to detect if the graph contains a negative cycle.

---

## Q2. Minimum Spanning Tree

Implement the Kruskal's and Prim's algorithms to find the minimum spanning tree of a graph.

### 2a. Prim's Algorithm

Implement the [Prim's algorithm](https://fahadsultan.com/csc223/datastructs/graphs_prim.html) to find the minimum spanning tree of a graph using a priority queue.

### 2b. Kruskal's Algorithm

Implement the [Kruskal's algorithm](https://fahadsultan.com/csc223/datastructs/graphs_kruskal.html) to find the minimum spanning tree of a graph. 

---

Feel free to define any additional <u>private</u> helper methods or classes as needed. 

You are to write JUnit test cases in corresponding test classes for each of the four classes above. Make sure to test all the methods keeping in mind any <u>edge cases</u> that may arise.

---

## **Instructions**
1. Click on `Accept Assignment` at the link posted on Moodle. This is equivalent to forking the repository.</font>
2. Navigate to the repository you just forked. It should be named `pa9-<your-github-username>`
3. Click on the green `Code` button. Click on SSH tab and copy the SSH link (it is important to use SSH and not HTTPS).
4. Open the terminal on your computer and navigate to the folder where you want to clone the repository.
5. Run `git clone <copied-ssh-link>`. Replace `<copied-ssh-link>` with the link you copied in step 3.
6. In Visual Studio Code, click on File -> Open and open the folder of the repository. 

## Evaluation

Your submission will be evaluated based on the correctness **_and completeness_** of your implementation and test cases. If your implementation is correct but your test cases are incomplete, you will lose points. If your implementation is incorrect, you will lose even more points.

## **Code Review**
This assignment will only be graded after you schedule and complete an in-person code review 
with the instructor within 5 working days after the due date, if not earlier. \
Use [this calendly link](https://calendly.com/ssultan-dpq/) provided on the course website to schedule a code review.

---
