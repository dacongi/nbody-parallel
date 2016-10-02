#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tree.h"
#include "body.h"
#include "mpi.h"

#define FILENAME "input_file_8.txt"
#define MASTER 0


//#define GENERATE_INPUT_FILE 1
#define DEBUG_GET_ARGUMENT 1
//#define READ_GENERATE_BODIES 1
//#define DEBUG_READ_FILE 1
#define DEBUG_OUTPUT_STATE 1
#define GENERATE_OUTPUT_FILE 1



//const double G = 6.673e-11;
const double G = 6;

const double treeratio = 1;


int my_rank;                            /* Rank del processo            */
int size;                               /* Numero dei processo          */



//prototipi di funzione
void Help();
void Get_Input_Arguments(int argc, char* argv[], int* num_particles, int* num_steps, double* delta_t, int* output_freq, char* init_condition);
struct body * Read_File_Init_Conditons(char *filename, struct body * bodies, int num_particles);
struct body * Generate_Init_Conditions(const int num_particles);
void Output_State(double time, struct body * bodies, int num_particles);
void Generate_Output_File(struct body *bodies, int num_particles);


double min (double a, double b);
double max (double a, double b);



//main program
int main(int argc, char *argv[]){
    
    
    int num_particles;                  /* Numero totale di particelle       */

    int steps;                          /* Step corrente                     */
    int num_steps;                      /* Numero timesteps                  */
    double delta_t;                     /* Taglia timestep                   */
    int output_freq;                    /* Frequenza di output               */
    char init_condition;                /* Condizione iniziale di input      */
    
    struct body *bodies;

    double time;                        /* Tempo corrente                    */
    double start_time;                  /* Tempo inizio esecuzione MPI       */
    double end_time;                    /* Tempo fine esecuzione MPI         */

    
    

    //INIZIALIZZO MPI
    MPI_Init(&argc, &argv);
    //CONSIDERO NUMERO DI PROCESSI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //CONSIDERO IL RANK DEL PROCESSO
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    //PRELEVO PARAMETRI DI INPUT DA LINEA DI COMANDO
    Get_Input_Arguments(argc,argv, &num_particles, &num_steps, &delta_t, &output_freq, &init_condition );
    
    
    //GENERO O LEGGO I VALORI INIZIALI DELLA SIMULAZIONE
    if(my_rank==MASTER){
        if (init_condition == 'g')
            bodies=Generate_Init_Conditions(num_particles);
        else
            bodies=Read_File_Init_Conditons(FILENAME,bodies,num_particles);
    
    }
    
    MPI_Bcast(&num_particles,1,MPI_INT, 0, MPI_COMM_WORLD);
    //printf("Sono il processo %d - num_particles: %d\n", my_rank,num_particles);
    
    if(my_rank!=MASTER){
        if(!(bodies=malloc(num_particles*sizeof(struct body)))){
            printf("Impossibile allocare memoria per %d bodies", num_particles);
            exit(1);
        }
    }
    
//    if(my_rank==MASTER)
//       printf("%lf %lf", bodies[2].m, bodies[1].vx);
    
    int i;
    for(i=0; i<num_particles; i++){
        MPI_Bcast(&((bodies)[i].m),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].x),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].y),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].fx),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].fy),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    
//    if(my_rank!=MASTER)
//        printf("%lf %lf", bodies[2].m, bodies[1].vx);
    

    
    start_time = MPI_Wtime();

    
#ifdef DEBUG_OUTPUT_STATE
    if(my_rank==MASTER)
        Output_State(time,bodies,num_particles);
#endif
    

    int step;
    //itero per il numero di step indicato
    for(step=1; step<=num_steps;step++){
        time = step*delta_t;
        
        struct node * rootnode;
        double xmin=0.0, xmax=0.0;
        double ymin=0.0, ymax=0.0;
        
        
        for(i=0; i<num_particles; i++){
            bodies[i].fx=0.0;
            bodies[i].fy=0.0;
            xmin=min(xmin,bodies[i].x);
            xmax=max(xmax,bodies[i].x);
            ymin=min(ymin,bodies[i].y);
            ymax=max(ymax,bodies[i].y);
        }
        
        rootnode=Create_Node(bodies+0,xmin,xmax,ymin,ymax);
        
        for(i=1; i<num_particles; i++)
            Insert_Body(bodies+i,rootnode);
        
        for(i=my_rank; i<num_particles; i+=size){
            Tree_Sum(rootnode,bodies+i, G,treeratio);
        }
        
        for(i=0; i<num_particles; i++){
            MPI_Bcast(&(bodies[i].fx),1,MPI_DOUBLE,i%size,MPI_COMM_WORLD);
            MPI_Bcast(&(bodies[i].fy),1,MPI_DOUBLE,i%size,MPI_COMM_WORLD);
        }
        
        //Aggiorno posizioni e velocità
        for(i=0; i<num_particles; i++){
            bodies[i].x += delta_t*bodies[i].vx;
            bodies[i].y += delta_t*bodies[i].vy;
            
            bodies[i].vx += bodies[i].fx * (delta_t/bodies[i].m);
            bodies[i].vy += bodies[i].fy * (delta_t/bodies[i].m);
            
        }

        
        Destroy_Tree(rootnode);
        
        
#ifdef DEBUG_OUTPUT_STATE
        if (steps % output_freq == 0)
            Output_State(time, bodies, num_particles);
#endif
        
    }
    

    
#ifdef GENERATE_OUTPUT_FILE
    if(my_rank==MASTER)
        Generate_Output_File(bodies, num_particles);
#endif
    
    
    end_time=MPI_Wtime();
    if(my_rank==MASTER){
        printf("Tempo trascorso = %e seconds\n", end_time-start_time);
        fprintf(stderr, "%d;%d;%e\n", num_particles,size,end_time-start_time );
    }
    
    
    
    

    
    free(bodies);
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
/* ./Help */

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
 *    init_condition:   da dove leggere input
 */
void Get_Input_Arguments(int argc, char* argv[], int* num_particles, int* num_steps, double* delta_t, int* output_freq, char* init_condition){
    

    if(argc!=6)
        Help();
        
    *num_particles = strtol(argv[1], NULL, 10);
    *num_steps = strtol(argv[2], NULL, 10);
    *delta_t = strtod(argv[3], NULL);
    *output_freq = strtol(argv[4], NULL, 10);
    *init_condition = argv[5][0];
    
    
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
    if (my_rank == MASTER) {
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
 * Out args:
 *    bodies:               array di strutture body
 */
struct body * Generate_Init_Conditions(const int num_particles)
{
    struct body * bodies;
    double mass = 5;
    double gap= 10;
    int i;
    
    //alloca memoria per contenere le particelle
    if(!(bodies = malloc(num_particles*sizeof(struct body))))
    {
        printf("Impossibile allocare memoria per %d particelle.\n",num_particles);
        return NULL;
    }
    
    for(i = 0; i < num_particles; i++)  //initialize with random position
    {
        bodies[i].m = mass;
        bodies[i].vx = 0.0;
        bodies[i].vy = 0.0;
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
        
        
        bodies[i].x = i*gap;
        bodies[i].y = 0.0;
        
#ifdef READ_GENERATE_BODIES
        printf("body %d: x=%f y=%f\n",i,bodies[i].x,bodies[i].y);
#endif
    }
    
#ifdef GENERATE_INPUT_FILE
    FILE *fp_generate = fopen("generate_input.txt","w+");
    if(!fp_generate){
        printf("Inpossibile creare il file <%s>\n", "generate_input");
        exit(1);
    }
    
    for(i=0; i<num_particles; i++){
        fprintf(fp_generate, "%lf %lf %lf %lf %lf\n", bodies[i].m, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy );
    }
    
    fclose(fp_generate);
#endif

    
    return bodies;
    
}
/* ./Generate_Init_Conditions */


/*---------------------------------------------------------------------
 * Function:   Read_File_Init_Conditons
 * Purpose:    Legge da un file le infomazioni inziali per ogni particella:
 massa, posizione x, posizione y, velocita x, velocita y
 * In args:
 *    filename:             nome del file in input
 *    num_particles:        numero totale di particelle
 * Out args:
 *    bodies:               array di strutture body
 */
struct body * Read_File_Init_Conditons(char *filename, struct body * bodies, int num_particles){
    int i;
    
    FILE * fp = fopen(filename, "r");
    if(!fp)
    {
        printf("Impossibile aprire il file <%s>\n\n",filename);
        exit(1);
    }
    
    
    if(!((bodies)=malloc(num_particles*sizeof(struct body))))
    {
        printf("Impossibile allocare memoria per %d bodies", num_particles);
        exit(1);
    }
    
    for(i = 0; i < num_particles; i++)
    {
        fscanf(fp,"%le %le %le %le %le\n", &((bodies)[i].m), &((bodies)[i].x), &((bodies)[i].y), &((bodies)[i].vx), &((bodies)[i].vy));
 
    }
    fflush(fp);
    fclose(fp);
    
    
    //printf("\n%d\n",num_particles);
    
#ifdef DEBUG_READ_FILE
    for(i=0; i<num_particles; i++){
        //printf("%d %d\n", i, num_particles);
        printf("%lf %lf %lf %lf %lf\n", bodies[i].m, bodies[i].x,bodies[i].y,bodies[i].vx,bodies[i].vy);
    }
#endif
    
    return bodies;
}
/* ./Read_File_Init_Conditons */


/*---------------------------------------------------------------------
 * Function:   Output_State
 * Purpose:    Stampa lo stato corrente del sistema
 * In args:
 *    time:                 current time
 *    bodies:               array di tutte le strutture body
 *    num_particles:        numero totale di particelle
 */
void Output_State(double time, struct body * bodies, int num_particles){
    int i;
    
    if (my_rank == MASTER) {
        printf("Current time:%.2f\n", time);
        for(i=0; i<num_particles; i++){
            printf("Particle:%d\tX:%f", i, bodies[i].x);
            printf("\tY:%f", bodies[i].y);
            printf("\tVx:%f ", bodies[i].vx);
            printf("\tVy:%f\n", bodies[i].vy);
        }
        printf("\n");
    }
}
/* ./Output_State*/


/*---------------------------------------------------------------------
 * Function:  Generate_Output_File
 * Purpose:   Genero un file di output con i risultati di ogni particella a termine della simulazione
 * In args:
 *    bodies:           array di tutte le strutture body
 *    num_particles:    numero totale di particelle
 */
void Generate_Output_File(struct body *bodies, int num_particles){
    int i;
    FILE *fp_generate_output = fopen("generate_output.txt","w+");
    if(!fp_generate_output){
        printf("Cannot create output file <%s>\n", "generate_output");
        exit(1);
    }
    
    for(i=0; i<num_particles;i++)
        fprintf(fp_generate_output, "%lf %lf %lf %lf %lf\n", bodies[i].m, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy );

    fflush(fp_generate_output);
    fclose(fp_generate_output);
}
/* ./Generate_Output_File*/
















double max (double a, double b)
{
    return (a > b ? a : b);
}

double min (double a, double b)
{
    return (a < b ? a : b);
}
