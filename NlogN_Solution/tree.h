#ifndef TREE
#define TREE

struct node {
    double totalmass;
    double centerx, centery;
    double xmin, xmax;
    double ymin, ymax;
    double diag;
    struct body * bodyp;
    struct node * NO;
    struct node * NW;
    struct node * SW;
    struct node * SO;
};




#endif




struct node *Create_Node(struct body * bodyp, double xmin, double xmax, double ymin, double ymax);
void Insert_Body(struct body * insbody, struct node * nodep);
void Tree_Sum(struct node * nodep, struct body * bodyp, double G, double ratiothreshold );
void Destroy_Tree(struct node * nodep);