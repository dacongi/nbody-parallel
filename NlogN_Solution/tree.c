#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "body.h"
#include "tree.h"

enum quadrant{NO,NW,SW,SO};


//firme
enum quadrant Get_Quadrant(double x, double y, double xmin, double xmax, double ymin, double ymax);
void Update_Center_Mass(struct node * nodep, struct body * bodyp);

















//CREO UN QUADRANTE (XMIN,XMAX,YMIN,YMAX) RITORNO IN QUALE SOTTOQUADRANTE IN CUI APPARTIENE (X,Y)
enum quadrant Get_Quadrant(double x, double y, double xmin, double xmax, double ymin, double ymax){
    
    double midx, midy;
    
    //dati i limiti del quadrante ne prendo i punt medi
    midx = xmin + 0.5*(xmax-xmin);
    midy = ymin + 0.5*(ymax-ymin);
    
    if(y>midy){
        if(x>midx)
            return NO;
        else
            return NW;
    }
    else{
        if(x>midx)
            return SO;
        else
            return SW;
    }

}

//CREO UN NODO FOGLIA DA INSERIRE ALL'INTERNO DELL'ALBERO(quadrante)
struct node *Create_Node(struct body *bodyp, double xmin, double xmax, double ymin, double ymax){
    struct node *rootnode;
    if(!(rootnode=malloc(sizeof(struct node)))){
        printf("Impossibile allocare il nodo\n");
    }
    
    rootnode->totalmass = bodyp->m;
    rootnode->centerx = bodyp->x;
    rootnode->centery = bodyp->y;
    rootnode->xmin = xmin;
    rootnode->xmax = xmax;
    rootnode->ymin = ymin;
    rootnode->ymax = ymax;
    
    rootnode->diag = sqrt(( pow(xmax - xmin, 2) + pow(ymax - ymin, 2) ));
    
    rootnode->bodyp = bodyp;
    
    //E' un nodo foglia nn si riferisce a nessun altro nodo;
    rootnode->NO = NULL;
    rootnode->NW = NULL;
    rootnode->SW = NULL;
    rootnode->SO = NULL;
        
    return rootnode;
}

//INSERISCO UNA PARTICELLA NELL'ALBERO, TRASFORMANDO UN NODO FOGLIA IN UN RAMO
void Insert_Body(struct body *insbody, struct node *nodep){
    
    enum quadrant existingquad, newquad;
    double xmid, ymid;
    
    xmid = nodep->xmin + 0.5*(nodep->xmax - nodep->xmin);
    ymid = nodep->ymin + 0.5*(nodep->ymax - nodep->ymin);
    
    
    //se il nodo è una foglia la trasformo in ramo ed inserisco la particella in un nuovo sotto quadrante
    if(nodep->bodyp != NULL){
        existingquad=Get_Quadrant(nodep->bodyp->x, nodep->bodyp->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
        
        switch (existingquad) {
            case NO:
                nodep->NO = Create_Node(nodep->bodyp, xmid, nodep->xmax, ymid, nodep->ymax);
                break;
            case NW:
                nodep->NW = Create_Node(nodep->bodyp, nodep->xmin, xmid, ymid, nodep->ymax);
                break;
            case SW:
                nodep->SW = Create_Node(nodep->bodyp, nodep->xmin, xmid, nodep->ymin, ymid);
                break;
            case SO:
                nodep->SO = Create_Node(nodep->bodyp, xmid, nodep->xmax, nodep->ymin, ymid);
                break;
        }
        
        nodep->bodyp = NULL;
    }

    
    newquad = Get_Quadrant(insbody->x, insbody->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
    Update_Center_Mass(nodep,insbody);
    
    //inserisco un nuovo punto in un nuovo quadrante se questo è vuoto altrimenti chiamo ricorsivamente Insert_Body
    switch (newquad){
        case NO:
            if(nodep->NO == NULL)
            {
                nodep->NO = Create_Node(insbody, xmid, nodep->xmax, ymid, nodep->ymax);
            } else {
                Insert_Body(insbody,nodep->NO);
            }
            break;
        case NW:
            if(nodep->NW == NULL)
            {
                nodep->NW = Create_Node(insbody, nodep->xmin, xmid, ymid, nodep->ymax);
            } else {
                Insert_Body(insbody,nodep->NW);
            }
            break;
        case SW:
            if(nodep->SW == NULL)
            {
                nodep->SW = Create_Node(insbody, nodep->xmin, xmid, nodep->ymin, ymid);
            } else {
                Insert_Body(insbody,nodep->SW);
            }
            break;
        case SO:
            if(nodep->SO == NULL)
            {			
                nodep->SO = Create_Node(insbody, xmid, nodep->xmax, nodep->ymin, ymid);
            } else {			
                Insert_Body(insbody,nodep->SO);
            }
            break;
    }
    
    
}



//AGGIORNO IL CENTRO DI MASSA IN SEGUITO AD UN INSERIMENTO
void Update_Center_Mass(struct node * nodep, struct body * bodyp)
{
    nodep->centerx = (nodep->totalmass*nodep->centerx + bodyp->m*bodyp->x)/(nodep->totalmass + bodyp->m);
    nodep->centery = (nodep->totalmass*nodep->centery + bodyp->m*bodyp->y)/(nodep->totalmass + bodyp->m);
    nodep->totalmass += bodyp->m;
    return;
}


//ELIMINO ALBERO RICORSIVAMENTE
void Destroy_Tree(struct node *nodep){
    if(nodep!=NULL){
        if(nodep->NO != NULL)
            Destroy_Tree(nodep->NO);
        if(nodep->NW != NULL)
            Destroy_Tree(nodep->NW);
        if(nodep->SW != NULL)
            Destroy_Tree(nodep->SW);
        if(nodep->SO != NULL)
            Destroy_Tree(nodep->SO);
        
        free(nodep);
    }
}


//SOMME LE FORZE SU BODYP
void Tree_Sum(struct node *nodep, struct body *bodyp, double G, double treeratio){
    double dx, dy, len, len_3;
    double m_g;
    double fact;
    
    dx = nodep->centerx - bodyp->x;
    dy = nodep->centery - bodyp->y;
    
    len = sqrt(pow(dx,2) + pow(dy,2));
    len_3 = pow(len,3);
    
    if( (((len / nodep->diag) > treeratio) || (nodep->bodyp))&&(nodep->bodyp!=bodyp) )
    {
        m_g = G*nodep->totalmass*bodyp->m;
        fact=m_g/len_3;
        
        
        
        bodyp->fx += fact*dx;
        bodyp->fy += fact*dy;
        
    } else {
        if(nodep->NO) { Tree_Sum(nodep->NO, bodyp, G, treeratio); }
        if(nodep->NW) { Tree_Sum(nodep->NW, bodyp, G, treeratio); }
        if(nodep->SW) { Tree_Sum(nodep->SW, bodyp, G, treeratio); }
        if(nodep->SO) { Tree_Sum(nodep->SO, bodyp, G, treeratio); }
    }
}




