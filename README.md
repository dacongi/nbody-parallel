# N-body simulation 
Examples MPI nbody problem implementation, project course of Parallel and Concurrent Programming.

### Problem statement

The problem is described [here](https://en.wikipedia.org/wiki/N-body_simulation).


### N^2 Solution

About the problem described in the “problem statement”, in this document’s section i described a basic parallel MPI based implementation of “n-body”.
Basic algorithm’s idea is that the only communication among the process/tasks occurs when the algorithm is computing the forces and, in order of the computing, each process/task needs the position and mass of every other particle. 
In this way, i use [MPI_Allgather](http://www.mpich.org/static/docs/v3.2/www3/MPI_Allgather.html) that is expressly designed for this situation; in fact this MPI’s function gathers data from all processes and distributes it to all processes.
An other important observation is that i don’t collect data about a single particle (e.g. mass, velocity, position) into a single datatype struct. However, if i use this datatype in MPI implementation, i’ll need to use a derived datatype in the call to MPI_Allgather, and the communications with this datatype tend to be slower than communications with basic MPI datatype. For this reason, i use individual arrays for the masses, positions, and velocities. 

About the implementation, i suppose that the array ‘positions’ can store the position of all system's particles. Further, i define 'vectorMPI' that is an MPI datatype that stores two contiguous doubles. Also i suppose that 'num_particles' (number of all particles into the system) is evenly divisible by size(number of process) and chunk=num_particles/size.

In this implementation, i made the following choice:
  - each process stores the entire global array of particles masses;
  - each process only uses a single n-element array for position;
  - each process uses a pointer 'my_positions' that refers to the start of its block of positions. Thus, on process 0   my_positions = positions; on process 1 my_position = position + chunk, and, so on.

So process 0 reads all the initial conditions into three n-element arrays. Since i’m storing all the masses on each process, i broadcast masses. Also, since each process will need the global array of positions for the first computation of forces, i broadcast positions. However, velocities are only used locally for the updates to positions and velocities, so we scatter velocieties.

### NLog(N) Solution 

add here a description...
http://arborjs.org/docs/barnes-hut

### Results

add chart and result ..

### Author

add info about you
