# N-body simulation 
Examples MPI nbody problem implementation, project course of Parallel and Concurrent Programming.
Professor: Vittorio Scarano

### Problem statement

In an n-body problem, we need to find the positions and velocities of a collection of interacting particles over a period of time. For example, an astrophysicist might want to know the positions and velocities of a collection of stars, while a chemist might want to know the positions and velocities of a collection of molecules or atoms.
An n-body solver is a program that finds the solution to an n-body problem by simulating the behavior of the particles. The input to the problem is the mass, position, and velocity of each particle at the start of the simulation, and the output is typically the position and velocity of each particle at a sequence of user-specified times, or simply the position and velocity of each particle at the end of a user-specified time period.

The problem is described [here](https://en.wikipedia.org/wiki/N-body_simulation).



### N^2 Solution

About the problem described in the “problem statement”, in this document’s section i described a basic parallel MPI based implementation of “n-body problem”.
Basic algorithm’s idea is that the only communication among the process/tasks occurs when the algorithm is computing the forces and, in order of the computing, each process/task needs the position and mass of every other system's particle. 
In this way, i use [MPI_Allgather](http://www.mpich.org/static/docs/v3.2/www3/MPI_Allgather.html) that is expressly designed for this situation; in fact this MPI function gathers data from all processes and distributes it to all processes.
An other important observation is that i don’t collect data about a single particle (e.g. mass, velocity, position) into a single datatype struct. However, if i use this datatype in MPI implementation, i’ll need to use a derived datatype in the call to MPI_Allgather, and the communications with this datatype tend to be slower than communications with basic MPI datatype. For this reason, i use individual arrays for the masses, positions, and velocities. [[1]](https://books.google.it/books?id=SEmfraJjvfwC&printsec=frontcover&hl=it&source=gbs_ge_summary_r&cad=0#v=onepage&q&f=false). Other importart reason why, in my implementation, i don't use struct datatype is that i haven't constraints the simulation in two-dimensions or three-dimensions. In fact using arrays i can simulate the interaction of the particles in space or in the plane only change a macro; if i use a struct i'd rewrite the entire code. 

About the implementation, i suppose that the array ‘positions’ can store the position of all system's particles. Further, i define 'vectorMPI' that is an MPI datatype that stores two contiguous doubles. Also i suppose that 'num_particles' (number of all particles into the system) is evenly divisible by size(number of process) and chunk=num_particles/size.

In this implementation, i made the following choice:
  - each process stores the entire global array of particles masses;
  - each process only uses a single n-element array for position;
  - each process uses a pointer 'my_positions' that refers to the start of its block of positions. Thus, on process 0   my_positions = positions; on process 1 my_position = position + chunk, and, so on.

It will also read the input and print the results.
So process 0 reads all the initial conditions into three n-element arrays. Since i’m storing all the masses on each process, i broadcast masses. Also, since each process will need the global array of positions for the first computation of forces, i broadcast positions. However, velocities are only used locally for the updates to positions and velocities, so we scatter velocieties.

This is a pseudocode for the basic solution:
[[1]](https://books.google.it/books?id=SEmfraJjvfwC&printsec=frontcover&hl=it&source=gbs_ge_summary_r&cad=0#v=onepage&q&f=false)
  - Get input data;
  - For each timestep:
    - if (timestep output) Print positions and velocities of particles,
    - for each local particle q -> Compute total force on q,
    - for each local particle q -> Compute position and velocity,
    - allgather local positions into global position array;
  - Print positions and velocities of all particles;  




### NLog(N) Solution 

About the problem described “problem statement”, in this document’s section i described  parallel MPI implementation of “n-body problem” using Barnes-Hut Tree-Code.
The Barnes-Hut algorithm, introduced by Josh Barnes and Piet Hut in 1986 [[2]](https://en.wikipedia.org/wiki/Barnes–Hut_simulation), describes a method to solve  “n-body problems”. The main difference with the previous implementation is that, instead of directly summing all the forces on a single particle, B-H’s implementation use a tree based approximation scheme to reduce the computational complexity from N^2 to NlogN.

About this algorithm, i highlight that the main idea is to group nearby particles and approximate them as a single particle. If the group is sufficiently far away, we can approximate its gravitational effects by using its center of mass [[3]](http://www.cs.princeton.edu/courses/archive/fall03/cs126/assignments/barnes-hut.html). The center of mass of a group of particles is the average position of a particle in that group, weighted by mass. For example, if two particles have position (x1,y1) and (x2,y2) and masses m1 and m2, then their total mass and center of mass (x, y) are given by: 

  m = m1 + m2;      
    x = (x1m1 + x2m2) / m;      
      y = (y1m1 + y2m2) / m.


The Barnes-Hut algorithm is, for these reasons, a scheme for grouping together particles that are sufficiently nearby. Analyzing the problem in a two-dimensional space, i divide, recursively, the set of particles into groups by storing them in a quad-tree. A quad-tree is similar to a binary tree, except that each node has 4 children (that in the implementation i described as NO, NW, SE, SO). Each node represents a region of the two dimensional space. 
In this way, the "root node" represents the whole space, and its four children represent the four main quadrants of the space. In fact, The tree is created by inserting the individual particles into tree through the main node. For this reason, when the particles are inserted into a tree node, tree node will divide the space represented by this node into four equal sized children nodes. Each child can in turn be broken into 4 subsquares to get its children, and so on. 

In this implementation, i made the following choice:
  - data struct node:
    - it represents a single node of the quad-tree;
  - data struct body:
    - it represents a single particle in the system.

This is a pseudocode for the basic solution: 
* Get input data;
* For each timestep:
    * if (timestep output) Print positions and velocities of particles,
    * tree construction;
    * force calculation
    * update particles;
* Print positions and velocities of all particles.

### Valitadion test
...
...

### Results

add chart and result ..

### Author

Marco Castaldo - 	Bachelor of Science degree in Computer Science - University of Salerno
