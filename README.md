# N-body simulation 
Examples MPI nbody problem implementation, project course of Parallel and Concurrent Programming.

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

1)This is a pseudocode for the basic solution:
2)  Get input data; 
  for each timestep {
    if (timestep output)
      Print positions and velocities of particles;
    for each local particle loc q 
      Compute total force on loc q;
    for each local particle loc q
      Compute position and velocity of loc q;
    Allgather local positions into global pos array;
  }
  Print positions and velocities of particles;



### NLog(N) Solution 

add here a description...
http://arborjs.org/docs/barnes-hut

### Results

add chart and result ..

### Author

add info about you
