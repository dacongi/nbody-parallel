#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define MASTER 0
#define DIM 2                           /* sistema bidimensionale       */
#define X 0                             /* coordinata x                 */
#define Y 1                             /* coordinata y                 */
#define FILENAME "input_file.txt"
#define DEBUG_GET_ARGUMENT 1
#define DEBUG_OUTPUT_STATE 1
//#define DEBUG_FORCES_BEFORE 0
//#define DEBUG_FORCES_AFTER 0
//#define DEBUG_READ_FILE 0
//#define DEBUG_UPDATE_BEFORE 0
//#define DEBUG_UPDATE_AFTER 0
//#define GENERATE_INPUT_FILE 1
//#define GENERATE_OUTPUT_FILE 1

typedef double vector[DIM];             /* Vettore di tipo double       */
//const double G = 6.673e-11;
const double G = 6;

int my_rank;                            /* Rank del processo            */
int size;                               /* Numero dei processo          */

MPI_Datatype vectorMPI;

vector *velocities = NULL;              /* Array velocità               */

//prototipi di funzione
void Help();
void Get_Input_Arguments(int argc, char* argv[], int* num_particles, int* num_steps, double* delta_t, int* output_freq, char* init_condition);
void Read_File_Init_Conditions(char* filaName, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void Generate_Init_Conditions(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void Output_State(double time, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void Compute_Force(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk);
void Update_Particles(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double delta_t);
void Generate_Output_File(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);

//main program
int main(int argc, char* argv[]){
    
    int num_particles;                  /* Numero totale di particelle       */
    int chunk;                          /* Numero particelle per rank        */
    int num_steps;                      /* Numero timesteps                  */
    int steps;                          /* Step corrente                     */
    int my_particles;                   /* Particella rank corrente          */
    int output_freq;                    /* Frequenza di output               */
    double delta_t;                     /* Taglia timestep                   */
    double time;                        /* Tempo corrente                    */
    double *masses;                     /* Array di tutte le masse           */
    vector* my_positions;               /* Array posizioni rank              */
    vector* positions;                  /* Array di tutte le posizioni       */
    vector* my_velocities;              /* Array velocità rank               */
    vector* my_forces;                  /* Array forze rank                  */
    
    char init_condition;                /* Condizione iniziale di input      */
    double start_time;                  /* Tempo inizio esecuzione MPI       */
    double end_time;                    /* Tempo fine esecuzione MPI         */
    
    //INIZIALIZZO MPI
    MPI_Init(&argc, &argv);
    //CONSIDERO NUMERO DI PROCESSI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //CONSIDERO IL RANK DEL PROCESSO
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //PRELEVO PARAMETRI DI INPUT DA LINEA DI COMANDO
    Get_Input_Arguments(argc, argv, &num_particles, &num_steps, &delta_t, &output_freq, &init_condition);
    //DIVIDO LE PARTICELLE PER IL NUMERO DI PROCESSI
    chunk=num_particles/size;
    
    masses = malloc(num_particles*sizeof(double));
    positions = malloc(num_particles*sizeof(vector));
    my_forces = malloc(chunk*sizeof(vector));
    my_positions = positions + my_rank*chunk;
    my_velocities = malloc(chunk*sizeof(vector));
    
    if (my_rank == MASTER)
        velocities = malloc(num_particles*sizeof(vector));
   
    //DEFINISCO UN TIPO DI DATIMPI CONTIGUO
    MPI_Type_contiguous(DIM, MPI_DOUBLE, &vectorMPI);
    MPI_Type_commit(&vectorMPI);
    
    //GENERO O LEGGO I VALORI INIZIALI DELLA SIMULAZIONE
    if (init_condition == 'f')
        Read_File_Init_Conditions(FILENAME, masses, positions, my_velocities, num_particles, chunk);
    else
        Generate_Init_Conditions(masses, positions, my_velocities, num_particles, chunk);

    
    start_time = MPI_Wtime();
    
    #ifdef DEBUG_OUTPUT_STATE
    Output_State(0.0, masses, positions, my_velocities, num_particles, chunk);
    #endif
    
    //ITERO PER IL NUMERO DI STEP INDICATO
    for (steps = 1; steps <= num_steps; steps++) {
        time = steps*delta_t;
        //COMPUTO LE FORZE PER OGNI PARTICELLA
        for (my_particles = 0; my_particles < chunk; my_particles++)
            Compute_Force(my_particles, masses, my_forces, positions, num_particles, chunk);
        //AGGIORNO I VALORI PER OGNI PARTICELLA
        for (my_particles = 0; my_particles < chunk; my_particles++)
            Update_Particles(my_particles, masses, my_forces, my_positions, my_velocities, num_particles, chunk, delta_t);
        
        
        MPI_Allgather(MPI_IN_PLACE, chunk, vectorMPI, positions, chunk, vectorMPI, MPI_COMM_WORLD);

        #ifdef DEBUG_OUTPUT_STATE
        if (steps % output_freq == 0)
            Output_State(time, masses, positions, my_velocities, num_particles, chunk);
        #endif
    }
    
    #ifdef GENERATE_OUTPUT_FILE
    Generate_Output_File(masses, positions,  my_velocities,  num_particles, chunk);
    #endif
    
    end_time=MPI_Wtime();
    if(my_rank==MASTER)
        printf("Tempo trascorso = %e seconds\n", end_time-start_time);

    MPI_Type_free(&vectorMPI);
    free(masses);
    free(positions);
    free(my_forces);
    free(my_velocities);
    if (my_rank == MASTER)
        free(velocities);
    
    //CHIUDO MPI
    MPI_Finalize();
    return 0;

}

//funzioni
/*---------------------------------------------------------------------
 * Function: Help
 * Purpose:  Stampa delle istruzioni per l'esecuzione del programma
 */
void Help(){
    fprintf(stderr, "   mpirun -np <numero processi> <nome programma>");
    fprintf(stderr, "   <numero particelle> <numbero di timesteps>");
    fprintf(stderr, "   <taglia timesteps> <output frequency>");
    fprintf(stderr, "   <f|g>");
    fprintf(stderr, "   'g': il programma genera le condizoni iniziali\n");
    fprintf(stderr, "   'f': leggo da file le condizioni iniziali\n");
    
    exit(0);
}
/* Help */

/*---------------------------------------------------------------------
 * Function:  Get_Input_Arguments
 * Purpose:   Prende i gli argomenti da linea di comando
 * In args:
 *    argc:            numero di argomenti
 *    argv:            array di argomenti
 * Out args:
 *    num_particles:    numero totale di particelle
 *    n_steps:          numero timesteps
 *    delta_t:          taglia timesteps
 *    output_freq:      frequenza output
 *    inir_condition:   da dove leggere input
 */
void Get_Input_Arguments(int argc, char* argv[], int* num_particles, int* num_steps, double* delta_t, int* output_freq, char* init_condition){
    
    if(my_rank==MASTER){
        if(argc!=6)
            Help();
        
        *num_particles = strtol(argv[1], NULL, 10);
        *num_steps = strtol(argv[2], NULL, 10);
        *delta_t = strtod(argv[3], NULL);
        *output_freq = strtol(argv[4], NULL, 10);
        *init_condition = argv[5][0];
    }
    
    MPI_Bcast(num_particles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(num_steps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(delta_t, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(output_freq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(init_condition, 1, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    
    
    if (num_particles <= 0 || *num_steps < 0 || *delta_t <= 0) {
        if (my_rank == MASTER)
            Help();
        MPI_Finalize();
        exit(0);
    }
    if (*init_condition != 'g' && *init_condition != 'f') {
        if (my_rank == MASTER)
            Help();
        MPI_Finalize();
        exit(0);
    }
    
#  ifdef DEBUG_GET_ARGUMENT
    if (my_rank == 0) {
        printf("num_particles = %d\n", *num_particles);
        printf("num_steps = %d\n", *num_steps);
        printf("delta_t = %e\n", *delta_t);
        printf("output_freq = %d\n", *output_freq);
        printf("init_conditon = %c\n", *init_condition);
    }
#  endif
}
/* ./Get_Input_Arguments */

/*---------------------------------------------------------------------
 * Function:  Generate_Init_Conditions
 * Purpose:   Genero le condizioni iniziali per ogni particella:
              massa, posizione x, posizione t, velocita x, velocita y
 * In args:
 *    num_particles:        numero totale di particelle
 *    chunk:                numero di particelle assegnate al processo
 * Out args:
 *    masses:               array di tutte le masse
 *    positions:            array di tutte le posizioni dell particelle
 *    my_velocites:         array di velocita del processo
 *
 * Note:      Genero le particella a distanza uguale sull'asse x con stessa massa e con velocità inziale pari a 0.
 *
 */
void Generate_Init_Conditions(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk){
    int part;
    int line;
    //double mass = 5.0e24;
    double mass = 5;
    //double gap = 1.0e5;
    double gap = 10;
    
    if (my_rank == MASTER) {
    
        for (part = 0; part < num_particles; part++) {
            masses[part] = mass;
            positions[part][X] = part*gap;
            positions[part][Y] = 0.0;
            velocities[part][X] = 0.0;
        }
        
        #ifdef GENERATE_INPUT_FILE
        FILE *fp_generate = fopen("generate_input","w+");
        if(!fp_generate){
            printf("Cannot create input file <%s>\n", "generate_input");
            exit(1);
        }
        
        for(line=0; line<num_particles; line++){
            fprintf(fp_generate, "%lf %lf %lf %lf %lf\n", masses[line], positions[line][X], positions[line][Y], velocities[line][X], velocities[line][Y]);
        }
        
        fclose(fp_generate);
        #endif

        
    }
    
    
    
    
    
    MPI_Bcast(masses, num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(positions, num_particles, vectorMPI, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(velocities, chunk, vectorMPI, my_velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
}
/* ./Generate_Init_Conditions*/

/*---------------------------------------------------------------------
 * Function:   Read_File_Init_Conditons
 * Purpose:    Legge da un file le infomazioni inziali per ogni particella:
               massa, posizione x, posizione y, velocita x, velocita y
 * In args:
 *    filename:             nome del file in input
 *    num_particles:        numero totale di particelle
 *    chunk:                numero di particelle assegnate al processo
 * Out args:
 *    masses:               array di tutte le masse
 *    positions:            array di tutte le posizioni dell particelle
 *    my_velocites:         array di velocita del processo
 */
void Read_File_Init_Conditions(char * fileName, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk){
    int line=0;
    
    if(my_rank == MASTER){
        FILE *fp = fopen(FILENAME,"r");
        if(!fp){
            printf("Cannot open input file <%s>\n", FILENAME);
            exit(1);
        }
        
        for(line=0; line<num_particles; line++){
            fscanf(fp, "%lf %lf %lf %lf %lf", &masses[line], &positions[line][X], &positions[line][Y], &velocities[line][X], &velocities[line][Y]);
        }
        
        #ifdef DEBUG_READ_FILE
        for(line=0; line<num_particles; line ++){
            printf("%lf %lf %lf %lf %lf\n", masses[line], positions[line][X], positions[line][Y], velocities[line][X], velocities[line][Y]);
        }
        #endif
        
        fclose(fp);
    }
    
    MPI_Bcast(masses, num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(positions, num_particles, vectorMPI, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(velocities, chunk, vectorMPI, my_velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
}
/* ./Read_File_Init_Conditions*/

/*---------------------------------------------------------------------
 * Function:   Output_State
 * Purpose:    Stampa lo stato corrente del sistema
 * In args:
 *    time:                 current time
 *    masses:               array di tutte le masse
 *    positions:            array di tutte le posizioni dell particelle
 *    my_velocites:         array di velocita del processo
 *    num_particles:        numero totale di particelle
 *    chunk:                numero di particelle assegnate al processo
 */
void Output_State(double time, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk) {
    int part;

    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
    if (my_rank == MASTER) {
        printf("Current time:%.2f\n", time);
        for (part = 0; part < num_particles; part++) {
            //       printf("%.3f ", masses[part]);
            printf("Particle:%3d\tX:%10.3e ", part, positions[part][X]);
            printf("\tY:%10.3e ", positions[part][Y]);
            printf("\tVx:%10.3e ", velocities[part][X]);
            printf("\tVy:%10.3e\n", velocities[part][Y]);
        }
        printf("\n");
    }
}
/* ./Output_state */

/*---------------------------------------------------------------------
 * Function:       Compute_Forces
 * Purpose:        Calcolo delle forze delle particelle in my_particles.
 * In args:
 *    my_particles:     particella corrente
 *    masses:           array di tutte le masse
 *    positions:        array di tutte le posizioni dell particelle
 *    num_particles:    numero totale di particelle
 *    chunk:            numero di particelle assegnate al processo
 * Out arg:
 *    my_forces:        array di forze del processo my particles
 */
void Compute_Force(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk){
    
    int k, part;
    double m_g;
    vector f_part_k;
    double len, len_3, fact;
    
    /* Indice corrispondente alla particelle locali */
    part = my_rank*chunk + my_particles;
    my_forces[my_particles][X] = my_forces[my_particles][Y] = 0.0;
    
    #ifdef DEBUG_FORCES_BEFORE
    printf("Proc %d > Current total force on part %d = (%.3e, %.3e)\n", my_rank, part, my_forces[my_particles][X], my_forces[my_particles][Y]);
    #endif
    
    for (k = 0; k < num_particles; k++) {
        if (k != part) {
            f_part_k[X] = positions[part][X] - positions[k][X];
            f_part_k[Y] = positions[part][Y] - positions[k][Y];
            
            len=sqrt(pow(f_part_k[X],2)+pow(f_part_k[Y],2));
            //len = sqrt(f_part_k[X]*f_part_k[X] + f_part_k[Y]*f_part_k[Y]);
            len_3=pow(len,3);
            //len_3 = len*len*len;
            
            m_g = G*masses[part]*masses[k];
            fact = m_g/len_3;
            
            f_part_k[X] *= fact;
            f_part_k[Y] *= fact;
            
            #ifdef DEBUG_FORCES_AFTER
            printf("Proc %d > Force on part %d due to part %d = (%.3e, %.3e)\n", my_rank, part, k, f_part_k[X], f_part_k[Y]);
            #endif
            
            /* Forza totale sulla particella */
            my_forces[my_particles][X] += f_part_k[X];
            my_forces[my_particles][Y] += f_part_k[Y];
        }
    }

}
/* ./Compute_Force*/

/*---------------------------------------------------------------------
 * Function:  Update_Particles
 * Purpose:   Aggiorna la velocita e la posizone di ogni particella in my_particles
 * In args:
 *    my_particles:     particella corrente
 *    masses:           array di tutte le masse
 *    my_forces:        array di forze del processo
 *    num_particles:    numero totale di particelle
 *    chunk:            numero di particelle assegnate al processo
 *    delta_t:          step size
 *
 * In/out args:
 *    my_positions:     array di posizione del processo
 *    my_velocities:    array di velocita del processo
 *
 * Note:  This version uses Euler's method to update both the velocity
 *    and the position.
 */
void Update_Particles(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double delta_t) {
    
    int part;
    double fact;
    
    part = my_rank*chunk + my_particles;
    fact = delta_t/masses[part];

    #ifdef DEBUG_UPDATE_BEFORE
    printf("   Proc %d > Before update of %d:\n", my_rank, part);
    printf("   Position  = (%.3e, %.3e)\n", my_positions[my_particles][X], my_positions[my_particles][Y]);
    printf("   Velocity  = (%.3e, %.3e)\n", my_positions[my_particles][X], my_positions[my_particles][Y]);
    printf("   Net force = (%.3e, %.3e)\n", my_forces[my_particles][X], my_forces[my_particles][Y]);
    #endif
    
    my_positions[my_particles][X] += delta_t * my_velocities[my_particles][X];
    my_positions[my_particles][Y] += delta_t * my_velocities[my_particles][Y];
    my_velocities[my_particles][X] += fact * my_forces[my_particles][X];
    my_velocities[my_particles][Y] += fact * my_forces[my_particles][Y];
    
    #ifdef DEBUG_UPDATE_AFTER
    printf("Proc %d > Position of %d = (%.3e, %.3e), Velocity = (%.3e,%.3e)\n", my_rank, part, my_positions[my_particles][X], my_positions[my_particles][Y],
           my_velocities[my_particles][X], my_velocities[my_particles][Y]);
    #endif
}
/* ./Update_Particles*/


/*---------------------------------------------------------------------
 * Function:  Generate_Output_File
 * Purpose:   Genero un file di output con i risultati di ogni particella a termine della simulazione
 * In args:
 *    masses:           array di tutte le masse
 *    num_particles:    numero totale di particelle
 *    chunk:            numero di particelle assegnate al processo
 *    positions:        array di tutte le posizioni
 *
 * In/out args:
 *    my_velocities:    array di velocita del processo
 *
 */
void Generate_Output_File(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk){
    int line;
    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
    if(my_rank==MASTER){
        FILE *fp_generate_output = fopen("generate_output","w+");
        if(!fp_generate_output){
            printf("Cannot create output file <%s>\n", "generate_output");
            exit(1);
        }
        
        for(line=0; line<num_particles; line++){
            fprintf(fp_generate_output, "%lf %lf %lf %lf %lf\n", masses[line], positions[line][X], positions[line][Y], velocities[line][X], velocities[line][Y]);
        }
        
        fclose(fp_generate_output);

    }
}
/* ./Generate_Output_File*/
