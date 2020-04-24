#ifndef DMAP2_GUARD
#define DMAP2_GUARD 

#include <math.h>
#include <stdio.h>




#define ELEMSIZE 5
#define BUFFERSIZE 10
#define NOTFOUND 999999

typedef float rtype;
typedef int itype;
typedef unsigned int uint;






struct BoundaryGrid{
    rtype lowerCorner[3];
    rtype upperCorner[3];
    rtype dX[3];
    rtype* data;
    uint inds[3];
};
typedef struct BoundaryGrid BoundaryGrid;


struct BoundaryArray{
    rtype lowerCorner[BUFFERSIZE*3];
    rtype upperCorner[BUFFERSIZE*3];
    BoundaryGrid* grid;
    uint nBounds;
};
typedef struct BoundaryArray BoundaryArray;

static BoundaryArray* bA_ptr = NULL;

static BoundaryArray** bA_list = NULL;


void prepare_N_arrays(uint N){
    bA_list = malloc(sizeof(*bA_list)*N);
}

void store_current_ptr_to(uint N){
    bA_list[N] = bA_ptr;
}

void switch_to_pointer_N(uint N){
    bA_ptr = bA_list[N];
}

uint get_bounds(){
    return bA_ptr->nBounds;
}

float get_smallest_distance(){
    return bA_ptr->grid[get_bounds()-1].dX[0];
}




BoundaryArray* init_boundary_array(uint nBounds,
                                   uint* lowerCornerInd,
                                   uint* upperCornerInd,
                                   uint* inds,
                                    rtype* lowerCorner,
                                    rtype* upperCorner){
    
    /*
    printf("%u \n",nBounds);
    printf("%u %u %u \n",lowerCornerInd[0],lowerCornerInd[1],lowerCornerInd[2]);
    printf("%u %u %u \n",upperCornerInd[0],upperCornerInd[1],upperCornerInd[2]);
    printf("%u %u %u \n",inds[0],inds[1],inds[2]);
    printf("%f %f %f \n",lowerCorner[0],lowerCorner[1],lowerCorner[2]);
    printf("%f %f %f \n",upperCorner[0],upperCorner[1],upperCorner[2]);
    */

    BoundaryArray* bA = malloc(sizeof(*bA));
    bA->nBounds = nBounds;
    bA->grid = malloc(nBounds*sizeof(BoundaryGrid));

    if(BUFFERSIZE<=nBounds){
        printf("Increase buffer size and compile again");
        return bA;    
    }

    for(uint i=0; i<nBounds; ++i){
        BoundaryGrid* bg = bA->grid+i;
        uint memSize = inds[i]*inds[i]*inds[i]*ELEMSIZE*sizeof(*(bg->data));
        bg->data = malloc(memSize);  
        memset(bg->data, 0, memSize );    
        if(i==0){
            for(uint j=0;j<3;++j){
                bA->lowerCorner[i*3+j] = lowerCorner[j];
                bA->upperCorner[i*3+j] = upperCorner[j];
                bg->lowerCorner[j] = lowerCorner[j];
                bg->upperCorner[j] = upperCorner[j];
            }
        }else{
            for(uint j=0;j<3;++j){
                rtype dX = (bA->upperCorner[(i-1)*3+j]-bA->lowerCorner[(i-1)*3+j])/inds[i];
                bA->lowerCorner[i*3+j] = bA->lowerCorner[(i-1)*3+j]+lowerCornerInd[i]*dX;
                bA->upperCorner[i*3+j] = bA->upperCorner[(i-1)*3+j]-upperCornerInd[i]*dX;
                bg->lowerCorner[j] = bA->lowerCorner[i*3+j];
                bg->upperCorner[j] = bA->upperCorner[i*3+j];     
            }
        }

        for(int j =0;j<3;++j){
            bg->inds[j] = inds[i];
            bg->dX[j] = (bg->upperCorner[j]-bg->lowerCorner[j])/bg->inds[j];  
            //printf("%f ",bg->upperCorner[j]);
        }
        //printf("\n ");

    }

    return bA;
}



void init_boundary_array_PY(uint nBounds,
                                   uint* lowerCornerInd,
                                   uint* upperCornerInd,
                                   uint* inds,
                                    rtype* lowerCorner,
                                    rtype* upperCorner){
    bA_ptr = init_boundary_array(nBounds,lowerCornerInd,
                                 upperCornerInd,inds,
                                 lowerCorner,upperCorner);
}


inline void convert_from_1D_to_3D(uint* newIndx, uint sIndx, uint* inds){
    newIndx[0] = (uint)(sIndx*1.0/(inds[1]*inds[2])); 
    newIndx[1] = (uint)((sIndx-newIndx[0]*(inds[1]*inds[2]))*1.0/inds[2]);
    newIndx[2] = (uint)(sIndx-newIndx[1]*inds[2]-newIndx[0]*inds[1]*inds[2]);
}


uint traverse_bArray_PY(uint* indStore, rtype* pos){

    uint bAIndx = indStore[0];
    uint sIndx = indStore[1]; 

    uint* inds = bA_ptr->grid[bAIndx].inds;    
    rtype* dX = bA_ptr->grid[bAIndx].dX;   
    if(sIndx >= inds[0]*inds[1]*inds[2]) return 0;  


    uint newIndx[3];
    convert_from_1D_to_3D(newIndx, sIndx, inds);

    for(uint i =0;i<3;++i){
        pos[i] = bA_ptr->grid[bAIndx].lowerCorner[i]+(dX[i]/2.0+dX[i]*newIndx[i]);
    }   

    return sIndx+1;
}


void add_data_bArray(rtype* tmp_data, uint* indStore){
    uint bAIndx = indStore[0];
    uint sIndx = indStore[1]; 

    for(int i=0;i<ELEMSIZE;++i){
        bA_ptr->grid[bAIndx].data[sIndx*ELEMSIZE+i] += tmp_data[i]; 
    }
}


void get_data_bArray(rtype* tmp_data, uint* indStore){
    uint bAIndx = indStore[0];
    uint sIndx = indStore[1]; 

    for(int i=0;i<ELEMSIZE;++i){
        tmp_data[i] = bA_ptr->grid[bAIndx].data[sIndx*ELEMSIZE+i]; 
        //printf("%f ",tmp_data[i]);
    }
    //printf("\n");
}




void save_to_file(char* fname){
    BoundaryArray* ptr = bA_ptr;
    FILE * file= fopen(fname, "wb");
    if (file != NULL) {
        fwrite(ptr, sizeof(*ptr), 1, file);
        fwrite(ptr->grid, sizeof(BoundaryGrid), ptr->nBounds, file);
        for(uint i=0;i<ptr->nBounds;++i){
            BoundaryGrid* bG = ptr->grid+i; 
            fwrite(bG->data, sizeof(*(bG->data)),
                    bG->inds[0]*bG->inds[1]*bG->inds[2]*ELEMSIZE,file);
            printf("%f %f %f %f %f \n",bG->data[0],bG->data[1],bG->data[2],bG->data[3],bG->data[4]);
        }
        fclose(file);
    }


}




void load_from_file(char* fname){
    BoundaryArray* bA = malloc(sizeof(*bA));
    FILE * file= fopen(fname, "rb");
    if (file != NULL) {
        fread(bA, sizeof(*bA), 1, file);
        //printf("%f %f %f \n",bA->lowerCorner[6],bA->lowerCorner[7],bA->lowerCorner[8]);

        bA->grid = malloc(sizeof(*(bA->grid))*bA->nBounds);
        fread(bA->grid,sizeof(*(bA->grid)),bA->nBounds,file);
        for(uint i=0;i<bA->nBounds;++i){
            BoundaryGrid* bG = bA->grid+i;
            int inds = bG->inds[0]*bG->inds[1]*bG->inds[2]*ELEMSIZE;
            printf("%i %i %i \n",bG->inds[0],bG->inds[1],bG->inds[2]);
            bG->data = malloc(sizeof(*(bG->data))*inds);
            fread(bG->data, sizeof(*(bG->data)),inds,file);
            printf("%f %f %f %f %f\n",bG->data[0],bG->data[1],bG->data[2],bG->data[3],bG->data[4]);
            printf("%i \n",inds);
        }
        fclose(file);
    }
    bA_ptr = bA;
};


/*

struct BoundaryArray{
    rtype lowerCorner[BUFFERSIZE*3];
    rtype upperCorner[BUFFERSIZE*3];
    BoundaryGrid* grid;
    uint nBounds;
};
*/




void add_data_to_cached(rtype value, rtype* pos, rtype* cached_data, uint* indStore){
    
    uint* bAIndx = indStore+0;
    uint* sIndx = indStore+1;
    uint use_cache = indStore[2];
    //printf("check: %u %u %u",indStore[0],indStore[1],indStore[2]);


    if(use_cache==1){
        use_cache = in_cached(pos,cached_data);
    }

    if(use_cache==0){
        //printf("check");
        find_indx(pos, bAIndx, sIndx, cached_data);
        indStore[2] = 1;
    }

   
    // printf("bAIndx: %u",*bAIndx);
    if(*bAIndx==NOTFOUND) return;
   
    bA_ptr->grid[*bAIndx].data[(*sIndx)*ELEMSIZE+4] += value;
    
}


int in_cached(rtype* pos, rtype* cached_data){
    for(uint i=0;i<3;++i){
        if(cached_data[i]>pos[i] ||
                    pos[i]>=cached_data[i]+cached_data[3+i]){
            return 0;
        }
    }
    return 1;
}





int in_cached_indx(rtype* pos, uint* indx){
    BoundaryArray* bA = bA_ptr;
    rtype cached_data[6];
    uint newIndx[3];
    uint sIndx = indx[1];
    uint* inds = bA->grid[indx[0]].inds;

    convert_from_1D_to_3D(newIndx, sIndx, inds);

    for(uint i=0; i<3; ++i){
        rtype dX = bA->grid[indx[0]].dX[i];
        rtype lowerCorner = bA->grid[indx[0]].lowerCorner[i];
        cached_data[3+i] = dX;
        cached_data[i] = lowerCorner+newIndx[i]*dX;        
    }
    return in_cached(pos,cached_data);
}



void find_indx(rtype* pos, uint* bAIndx, uint* sIndx, rtype* cached_data){

    BoundaryArray* bA = bA_ptr;
    for(int i = bA->nBounds-1; i >= 0; --i){
        //printf("iter, %f %f %f %u\n",pos[0],pos[1],pos[2],i);
        rtype* lowerCorner = bA->lowerCorner + i*3;
        rtype* upperCorner = bA->upperCorner + i*3;
        
        int inside = 0;        
        for(uint j=0; j<3; ++j){
            if(lowerCorner[j]<=pos[j] && pos[j]<upperCorner[j]){
                inside++;
            }
            
        }
        if(inside==3){
            //printf("went in");
            *bAIndx = i;
            BoundaryGrid* bg = bA->grid+i;


            
            uint inds[3];
            for(uint j=0;j<3;++j){
                inds[j] = (uint)((pos[j]-lowerCorner[j])/bg->dX[j]);
                if(inds[j]>=bg->inds[j]){
                    inds[j] = bg->inds[j]-1;
                }     
            }


            for(uint k = 0; k<3; ++k){
                cached_data[0+k] = lowerCorner[0+k]+inds[k]*bg->dX[0+k];
                cached_data[3+k] = bg->dX[0+k];
            }


            *sIndx = (inds[0]*(bg->inds[1]*bg->inds[2]) + inds[1]*(bg->inds[2])+inds[2]);
            return;
        }
        //printf("haaa? \n");
    }
    //printf("not found \n");
    *bAIndx = NOTFOUND;
    *sIndx = NOTFOUND;
    //printf("not found \n");
}




void get_data_at(rtype* pos, rtype* data, rtype* cached_data, uint* indStore){
    uint use_cache = indStore[2];
    uint* bAIndx = indStore+0;
    uint* sIndx = indStore+1; 

    if(use_cache==1){
        use_cache = in_cached(pos,cached_data);
    }

    if(use_cache==0){;
        find_indx(pos, bAIndx, sIndx, cached_data);
        indStore[2] = 1;
    }else{
        return;
    }

    
    if(*bAIndx==NOTFOUND){
        //printf("not found");
        return;
    }
    for(uint i = 0; i<ELEMSIZE-1; ++i){
        data[i] = bA_ptr->grid[*bAIndx].data[(*sIndx)*ELEMSIZE+i];
    }  
    rtype dX =  bA_ptr->grid[*bAIndx].dX[0];
    rtype normalization = dX*dX*dX;
    data[ELEMSIZE-1] = bA_ptr->grid[*bAIndx].data[(*sIndx)*ELEMSIZE+ELEMSIZE-1]/normalization; 
    //printf("indstore: %u %u %u \n",indStore[0],indStore[1],indStore[2]); 
}





void get_data_at0(rtype* pos, rtype* data){
    BoundaryArray* bA = bA_ptr;
    for(int i = bA->nBounds-1; i >= 0; --i){
        rtype* lowerCorner = bA->lowerCorner + i*3;
        rtype* upperCorner = bA->upperCorner + i*3;
        int inside = 0;        
        for(uint j=0; j<3; ++j){
            if(lowerCorner[j]<pos[j] && upperCorner[j]>pos[j]){
                inside++;
            }
        }
        if(inside==3){
            BoundaryGrid* bg = bA->grid+i;
            uint inds[3];
            for(uint j=0;j<3;++j){
                inds[j] = (uint)((pos[j]-lowerCorner[j])/bg->dX[j]);            
            }
            uint ind = (inds[0]*(bg->inds[1]*bg->inds[2]) + inds[1]*(bg->inds[2])+inds[2])*ELEMSIZE;
            for(uint j=0;j<ELEMSIZE;++j){
                data[j] = bg->data[ind+j];
            }
            return;
        }
    }

}

#endif
