/*

This program was written by Miguel Arevalillo-Herráez and Hussein Abdul-Rahman 
to program the algorithm entitled the robust three-dimensional best path phase unwrapping algorithm that avoids singularity loops 

If you have any questions regarding this program or using the algorithm, please email 
miguel.arevalillo@uv.es, Hussein_Nemer@yahoo.com or m.a.gdeisat@ljmu.ac.uk

This program was written on 22nd April 2009

The wrapped phase volume is assumed to be of floating point data type. The resultant unwrapped phase volume is also of floating point type.
read the data from the file frame by frame.


To run the program do the following

	1- compile the C program. This will produce an executable file. Suppose that the name of this file is 3DBPPUASL.exe 
	(three-dimensional best path phase unwrapping avoiding singularity loops)

	2- Go to the DOS command prompot and type the following
		
3DBPPUASL [-t] wrapped_phase_input_file  unwrapped_phase_output_file  No_of_rows  No_of_columns  No_of_frames  [-m mask_map]  [-q] quality map

		wrapped_phase_input_file: The name of the file that contains the wrapped phase volume data. 
		The data are binary floating point type.
		
		unwrapped_phase_input_file: The name of the file that contains the unwrapped phase volume data. 
		The data produced by this program are binary floating point type.
		
		No_of_rows: The number of rows in the wrapped phase volume

		No_of_columns: The number of columns in the wrapped phase volume

		No_of_frames: The number of frames in the wrapped phase volume

		[] This means optional input

		-t use this switch if the input wrapped phase data type is ASCII format

		-m mask_map: Use this switch if you would like the algorithm to use the mask map supported to the algorithm by the user. 
		The mask map must have the same number of voxels (pixels) as the wrapped phase volume and it has an unsigned char data type format. 
		A value of 1 in the mask map corresponds to a valid voxel. Conversely, a value of 0 in the mask map corresponds to an invalid voxel.
		
		-q quality map: Use this switch if you would like to use your own quality map. 
		If you do not set this switch, the algorithm will use its own quality map. 
		The quality map must have the same number of voxels as the wrapped phase volume and it has floating point data format.
		Larger values in the quality map correspond to high quality voxels in the wrapped phase volume and vice versa. 

*/

/* includes */
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <float.h>

/* Constantes */
const float pi = 3.141592F;
const float twopi = pi*2;
const float MIN_RELIABILITY=1E4F;
int rows;
int columns;
int depth;
int number_of_pixels;
int twodsize;
int r_twodsize;
int debug=1;


/* pre-definiciones */
struct frontier;
struct residue_node;

/*************************************************************************************/
/*                                    ESTRUCTURAS                                    */  
/*************************************************************************************/

/* Estructura group */

// The structure to store each pixel in the image. It allows to store 
// information about the groups using a data structure which is
// similar to a linked list. The head of the group stores additional 
// information, such as its size and the number of borders with
// other groups it has.
struct group
{
  float value;               // The value of the pixel
  long current_increment;    // The number of 2pi increments with respect to the 
                             // pixels pointed by previous_group
  group* last_pixel;         // The last pixel of the group
  group* next_pixel;         // The next pixel of the group
  group* previous_group;     // A pixel in the same group
  long group_size;
  float reliability; // The reliability of the pixel
  int group_number;  // used to control the group the pixel belonga to
};

/* Estructura Distances */

// Temporal structure which is used in the calculation of the pixel reliability.
// The only purpose of the structure is to avoid repetitions in the calculations.
struct distances
{
   float horizontal_distance;
   float vertical_distance;
   float in_depth_distance;
   float diagonal_right_distance;
   float diagonal_left_distance;
   float diagonal_right_above_distance;
   float diagonal_left_above_distance;
   float diagonal_right_side_distance;
   float diagonal_left_side_distance;
   float opposite_corner_from_me;
   float opposite_corner_from_right;
   float opposite_corner_from_bottom;
   float opposite_corner_from_opposite;
   int horizontal_steps;
   int vertical_steps;
   int in_depth_steps;
};

// This structure is used to represent the residues. We use a 3d view of boxes. 
// which are represented by their centres. In this structure we store the sum of
// each connected four voxel group. Those which do not add to zero are the residues.
struct residue_values {
   float residue_front_sum;
   float residue_side_sum;
   float residue_top_sum;
   int treated[3]; // Indicates if the residues are already a part of a loop
};      

// This structure is used to store the paths. The path is considered as a line
// connecting the centers of the boxes
struct residue_node {
   int x;
   int y;
   int z;
   residue_values* value;
   residue_node* next;
   residue_node* back;
};

struct center_of_gravity_struct {
   float x;
   float y;
   float z;
   int number_of_elements;
};
       
// Structure that represents a border. It also stores the index in the heap the border
// is stored in. The number of 2pi jumps which would be required to locally unwrap
// the groups at both sides of the border is stored in steps
struct frontier
{
        group* group1;
        group* group2;
        float distance;
        int steps;
        float distance_backup;
        int is_banned;
};



/*************************************************************************************/
/*                                 GLOBAL VARIABLES                                  */  
/*************************************************************************************/

group *array;           // Stores the image and the groups
residue_values *residues;      // Stores the value of every four pixel loop
distances *derivatives; // It stores the first horizontal and vertical derivative. It is 
                        // temporarily used in the calculation of the pixel reliability
                        
frontier* frontier_array; // The array to store the borders
frontier* hfrontiers;     // Where the frontiers in x start
frontier* dfrontiers;     // where the frontiers in z start
unsigned char* mask;
float* qualities;

// the frontiers in y start in frontier_array
float* values_array; // To store the values load from the file temporally
int number_of_frontiers; // total number of frontiers

/*************************************************************************************/
/*************************************************************************************/
/*						                                                             */
/*                              PROCEDURES AND FUNCTIONS                             */  
/*						                                                             */
/*************************************************************************************/
/*************************************************************************************/

/*************************************************************************************/
/*                            FILE MANIPULATION FUNCTIONS                            */  
/*************************************************************************************/
// Loads a plain text file with the filename specified by the string 'name' in the global
// variable array, and stores the values in the array specified by temp. Allocation of temp
// must be done outside the function
void load_file_plain_text(char* name, float* temp)
{
   FILE *f;
   f=fopen(name,"r");
   //float* temp=values_array; 
   if (f==NULL)
     printf("ERROR - The file could not be opened in text mode");
   else
   {
     for (int i=0;i<number_of_pixels;i++){
          fscanf(f,"%f",temp);
          temp++;
     };
     fclose(f);
   }
}

// Loads a binary file with the filename specified by the string 'name' in the global
// variable array, and stores the values in the array specified by temp. Allocation of temp
// must be done outside the function
void load_file_binary(char* name, float* array)
{
   FILE *f;
   float temp;
	
   f=fopen(name,"rb");
   if (f==NULL) {
     printf("ERROR - Binary file could not be opened");
     exit(8);
   } else {
       fread(array,number_of_pixels*sizeof(float),1,f);
       fclose(f);
   }
}


// Loads a mask (binary) file with the filename specified by the string 'name' in the global
// variable array. The values are stored in the array specified as the second parameter. 
// The space allocation has to be performed in the main function.
void load_mask_plain_text(char* name, unsigned char* temp)
{
   FILE *f;
   f=fopen(name,"r");
   //float* temp=values_array; 
   if (f==NULL)
     printf("ERROR - The file could not be opened in text mode");
   else
   {
     for (int i=0;i<number_of_pixels;i++){
          fscanf(f,"%f",temp);
          temp++;
     };
     fclose(f);
   }
}
// Loads a mask (binary) file with the filename specified by the string 'name' in the global
// variable array. The values are stored in the array specified as the second parameter. 
// The space allocation has to be performed in the main function.
void load_mask_binary(char* name, unsigned char* array)
{
   FILE *f;
   float temp;
	
   f=fopen(name,"rb");
   if (f==NULL) {
     printf("ERROR - Binary file could not be opened");
     exit(8);
   } else {
       fread(array,number_of_pixels*sizeof(unsigned char),1,f);
       fclose(f);
   }
}


// Stores the signal in the memory in a file with the name specified by the string 'name'
// The file is stored in text format
void save_file_plain_text(char *name,group* array)
{
     group* p;
     p=array;
     float value;
     FILE* f;
     f=fopen(name,"w");
     for (int i=0;i<number_of_pixels;i++)
     {
       value=(p->value)+(p->current_increment)*twopi;
       fprintf(f,"%f\n",value);
       p++;
     }
     fclose(f);
}

// Stores the signal in the memory in a file with the name specified by the string 'name'
// The file is stores in binary format
void save_file_binary(char *name,group* array)
{
     group* p;
     p=array;
     float value;
     FILE* f;
     f=fopen(name,"wb");
     for (int i=0;i<number_of_pixels;i++)
     {
       value=(p->value)+(p->current_increment)*twopi;
       fwrite(&value,sizeof(float),1,f);
       p++;
     }
     fclose(f);
}

// This function is used to tranfer the float data in values_array (the values
// of the original signal) into the array of structs which is used by the algorithm
void convert_values_to_struct(){
   group* p=array;
   FILE *f;
   float* temp=values_array; 
   for (int i=0;i<number_of_pixels;i++) {
     p->value=*temp;
     p->next_pixel=NULL;
     p->last_pixel=p;
     p->current_increment=0;
     p->group_size=1;
     p++;
     temp++;
   };
}

// Create a new node from its pointer (the pointer determines the position)
// It takes a residue node and it calculates its coordinates, storing them in
// a residue_node struct and returning a pointer to this structure
residue_node* fill_node(residue_values* residue) {
  residue_node* this_node;
  this_node=new residue_node;
  int temp=residue-residues;
  this_node->z=temp/twodsize;
  temp=temp%twodsize;
  this_node->y=temp/columns;
  this_node->x=temp%columns;
  this_node->next=NULL;
  this_node->back=NULL;
  this_node->value=residue;
  return this_node;
}
      
// It follows the residues which enter the node to build a path 
// The path it builds is a double linked list
// It follows the path until either it finds a border or it
// returns to the starting point
// It the path does not end at a border, it builds the entire path. Otherwise,
// it builds the section of the path from the residue backward until it reaches 
// the border.
// When a residue is annotated as part of a loop, its treated flag is set to one
// to avoid that it is taken more than once. 
// It returns a pointer to the last element
residue_node* build_loop_backward(residue_values* residue,residue_node* last_node) {
  residue_node* this_node;
  residue_node* new_node;

  if (residue->residue_front_sum==1 && residue->treated[0]==0) {
    residue->treated[0]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
    
    if (this_node->z!=0) 
       return build_loop_backward(residue-twodsize,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue-twodsize;
       new_node->z=-1;
       new_node->next=this_node;
       return new_node;
    }
  }

  if (residue->residue_side_sum==1 && residue->treated[1]==0) {
    residue->treated[1]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
  
    if (this_node->x!=0) 
       return build_loop_backward(residue-1,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue-1;
       new_node->x=-1;
       new_node->next=this_node;
       return new_node;
    }   
   }

  if (residue->residue_top_sum==1 && residue->treated[2]==0) {
    residue->treated[2]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
  
    if (this_node->y!=0) 
        return build_loop_backward(residue-columns,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue-columns;
       new_node->y=-1;
       new_node->next=this_node;
       return new_node;
    }
    this_node->next=last_node;
    return this_node;
   }

  if ((residue+1)->residue_side_sum==-1 && (residue+1)->treated[1]==0) {
    (residue+1)->treated[1]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
  
    if (this_node->x<columns-2) 
       return build_loop_backward(residue+1,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue+1;
       new_node->x=columns-1;
       new_node->next=this_node;
       return new_node;
    }
    this_node->next=last_node;
    return this_node;
   }
  
   if ((residue+twodsize)->residue_front_sum==-1 && (residue+twodsize)->treated[0]==0) {
    (residue+twodsize)->treated[0]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
  
    if (this_node->z<depth-2) 
       return build_loop_backward(residue+twodsize,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue+twodsize;
       new_node->z=depth-1;
       new_node->next=this_node;
       return new_node;
    }

    this_node->next=last_node;
    return this_node;
   }

   if ((residue+columns)->residue_top_sum==-1 && (residue+columns)->treated[2]==0) {
    (residue+columns)->treated[2]=1; // We annotate that it has already been 
                                    // treated
    this_node=fill_node(residue);  // Builds a new residue node
    this_node->next=last_node;
    if (last_node!=NULL)
      last_node->back=this_node;
  
    if (this_node->y<rows-2) 
       return build_loop_backward(residue+columns,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->back=new_node;
       new_node->value=residue+columns;
       new_node->y=rows-1;
       new_node->next=this_node;
       return new_node;
    }
    this_node->next=last_node;
    return this_node;
   }  
 
   return last_node;
 }                         



// It follows the residues which exit the node to build a path 
// The path it builds is a double linked list
// It follows the path until either it finds a border or it
// returns to the starting point
// It the path does not end at a border, it builds the entire path. Otherwise,
// it builds the section of the path from the residue forward until it reaches 
// the border
// When a residue is annotated as part of a loop, its treated flag is set to one
// to avoid that it is taken more than once. 


residue_node* build_loop_forward(residue_values* residue,residue_node* last_node) {
  residue_node* this_node;
  residue_node* new_node;

  if (residue->residue_front_sum==-1 && residue->treated[0]==0) {
    // The path continues towards the front of the cube                                   
    residue->treated[0]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->z!=0)  {
       this_node->next=build_loop_forward(residue-twodsize,this_node); // General case
    }  
       
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue-twodsize;
       new_node->z=-1;
       new_node->back=this_node;
    }
    this_node->back=last_node;
    return this_node;
   }

  if (residue->residue_side_sum==-1 && residue->treated[1]==0) {
    // The path continues towards the left of the cube                                   
    residue->treated[1]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->x!=0) 
       this_node->next=build_loop_forward(residue-1,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue-1;
       new_node->x=-1;
       new_node->back=this_node;
    }
    this_node->back=last_node;
    return this_node;
   }

  if (residue->residue_top_sum==-1 && residue->treated[2]==0) {
    // The path continues towards the top of the cube                                   
    residue->treated[2]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->y!=0) 
       this_node->next=build_loop_forward(residue-columns,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue-columns;
       new_node->y=-1;
       new_node->back=this_node;
    }
    this_node->back=last_node;
    return this_node;
   }

  if ((residue+1)->residue_side_sum==1 && (residue+1)->treated[1]==0) {
    // The path continues towards the right of the cube                                   
    (residue+1)->treated[1]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->x<columns-2) 
       this_node->next=build_loop_forward(residue+1,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue+1;
       new_node->x=columns-1;
       new_node->back=this_node;
    }
    this_node->back=last_node;
    return this_node;
   }
  
   if ((residue+twodsize)->residue_front_sum==1 && (residue+twodsize)->treated[0]==0) {
    // The path continues towards the back of the cube                                  
    (residue+twodsize)->treated[0]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->z<depth-2) 
       this_node->next=build_loop_forward(residue+twodsize,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue+twodsize;
       new_node->z=depth-1;
       new_node->back=this_node;
    }
    this_node->back=last_node;
    return this_node;
   }

   if ((residue+columns)->residue_top_sum==1 && (residue+columns)->treated[2]==0) {
    // The path continues towards the bottom of the cube                                  
    (residue+columns)->treated[2]=1; // We annotate that it has already been 
                           // treated
    this_node=fill_node(residue);  // Builds a residue node
  
    if (this_node->y<rows-2) 
       this_node->next=build_loop_forward(residue+columns,this_node); // General case
    else {
       new_node=fill_node(residue); 
       this_node->next=new_node;
       new_node->value=residue+columns;
       new_node->y=rows-1;
       new_node->back=this_node;
    }

    this_node->back=last_node;
    return this_node;

   }  

   return NULL; // There are no more residues
}

residue_node* build_loop_forward(residue_values* residue) {
              return build_loop_forward(residue,NULL);
} 

residue_node* build_loop_backward(residue_values* residue) {
              return build_loop_backward(residue,NULL);
} 



// Adds nodes in the X axis until valor is reached.
// int arrive -> one if valor is included, cero otherwise
residue_node* advance_in_x(residue_node* node,int valor,int arrive) {
   // node-> the node at which we start
   // int valor -> stop value
   // int arrive -> one if valor is included, cero otherwise
   int new_value;
   residue_node* new_node=node; // In case it does not enter the loop
   new_value=(node->x)+1;
   while (new_value<valor+arrive) { // We include the last element
       new_node=new residue_node;
       node->next=new_node;
       new_node->back=node;
       new_node->x=new_value;
       new_node->y=node->y;
       new_node->z=node->z;
       new_node->value=node->value+1;       
       new_value++;
       node=new_node;
     }
   new_value=(node->x)-1;
   while (new_value>valor-arrive) { // We include the last element
     new_node=new residue_node;
     node->next=new_node;
     new_node->back=node;
     new_node->x=new_value;
     new_node->y=node->y;
     new_node->z=node->z;
     new_node->value=node->value-1;       
     new_value--;
     node=new_node;
   }
   
   node->next=NULL;
   return new_node;
}
     
// Adds nodes in the Y axis until valor is reached.
// int arrive -> one if valor is included, cero otherwise
residue_node* advance_in_y(residue_node* node,int valor,int arrive) {
   // node-> the node at which we start
   // int valor -> stop value
   residue_node* new_node=node; // In case it does not enter the loop
   int new_value=(node->y)+1;
   while (new_value<valor+arrive) { 
     new_node=new residue_node;
     node->next=new_node;
     new_node->back=node;
     new_node->x=node->x;
     new_node->y=new_value;
     new_node->z=node->z;
     new_node->value=node->value+columns;       
     new_value++;
     node=new_node;
   }
   new_value=(node->y)-1;
   while (new_value>valor-arrive) { 
     new_node=new residue_node;
     node->next=new_node;
     new_node->back=node;
     new_node->x=node->x;
     new_node->y=new_value;
     new_node->z=node->z;
     new_node->value=node->value-columns;       
     new_value--;
     node=new_node;
   }
   node->next=NULL;
   return new_node;
}

// Adds nodes in the Z axis until valor is reached.
// int arrive -> one if valor is included, cero otherwise
residue_node* advance_in_z(residue_node* node,int valor,int arrive) {
   // node-> the node at which we start
   // int valor -> stop value
   residue_node* new_node=node; // In case it does not enter the loop
   int new_value=(node->z)+1;
   while (new_value<valor+arrive) { 
     new_node=new residue_node;
     node->next=new_node;
     new_node->back=node;
     new_node->x=node->x;
     new_node->y=node->y;
     new_node->z=new_value;
     new_node->value=node->value+twodsize;       
     new_value++;
     node=new_node;
   }
   new_value=(node->z)-1;
   while (new_value>valor-arrive) { 
     new_node=new residue_node;
     node->next=new_node;
     new_node->back=node;
     new_node->x=node->x;
     new_node->y=node->y;
     new_node->z=new_value;
     new_node->value=node->value-twodsize;       
     new_value--;
     node=new_node;
   }
   node->next=NULL;
   return new_node;
}      


// It completes a loop when the head and tail are not the same
// This is necessary when the path ends at the borders. Then both
// ends need to be joint
// It does nothing if the path is actually a loop
void complete_loop(residue_node* head,residue_node* tail) {
   int minimum_option;
   int minimum_value;
   int this_value;

          if ((head->x==-1) || (head->x==columns-1)) { // The head is at the left or right side
             if (tail->x==head->x) { // It is on the same side
                 head=advance_in_y(head,tail->y,0);
                 if (tail->y==head->y)
                   head=advance_in_z(head,(tail->z),0);
                 else 
                   head=advance_in_z(head,(tail->z),1);
             } else
                if (tail->y==-1) { // The tail is above
                  head=advance_in_y(head,-1,1);
                  head=advance_in_x(head,tail->x,0);
                 if (tail->x==head->x)
                    head=advance_in_z(head,tail->z,0);
                  else
                    head=advance_in_z(head,tail->z,1);
                } else
                    if (tail->y==rows-1) { // The tail is below
                       head=advance_in_y(head,rows-1,1);
                       head=advance_in_x(head,tail->x,0);
                       if (tail->x==head->x)
                         head=advance_in_z(head,tail->z,0);
                       else
                         head=advance_in_z(head,tail->z,1);
                    } else
                      if (tail->z==-1) { // The tail is in front 
                         head=advance_in_z(head,-1,1);
                         head=advance_in_x(head,(tail->x),0);
                         if (tail->x==head->x)
                            head=advance_in_y(head,tail->y,0);
                         else
                            head=advance_in_y(head,tail->y,1);
                      } else 
                         if (tail->z==depth-1) {
                          // The tail is at the back face
                         head=advance_in_z(head,depth-1,1);
                         head=advance_in_x(head,(tail->x),0);
                         if (tail->x==head->x)
                            head=advance_in_y(head,tail->y,0);
                         else
                            head=advance_in_y(head,tail->y,1);
                      } else {
                          // It is on the opposite side
                          minimum_option=0;
                          minimum_value=2*(depth-2)-(head->z)-(tail->z); //towards the end

                          this_value=(head->z)+(tail->z); // towards me
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=1;
                          }
             
                          this_value=(head->y)+(tail->y); // up
                          if (this_value<minimum_value){
                            minimum_value=this_value;
                            minimum_option=2;
                          }
             
                          this_value=2*(rows-2)-(head->y)-(tail->y); // down
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=3;
                          }
             
                          switch (minimum_option){
                              case 0: head=advance_in_z(head,depth-1,1);break;
                              case 1: head=advance_in_z(head,-1,1); break;
                              case 2: head=advance_in_y(head,-1,1); break;
                              case 3: head=advance_in_y(head,rows-1,1); break;
                         }
                         if (head->x==columns-1)
                              head=advance_in_x(head,-1,1);
                         else
                              head=advance_in_x(head,columns-1,1);
                         complete_loop(head,tail);
                      }
           } else
           if ((head->y==-1) || (head->y==rows-1)) { // The head is at the top or botton
             if (tail->y==head->y) { // It is on the same side
                 head=advance_in_x(head,tail->x,0);
                 if (tail->x==head->x)
                     head=advance_in_z(head,tail->z,0);
                 else
                     head=advance_in_z(head,tail->z,1);

             } else
                if (tail->x==-1) { // The tail is on the left
                  head=advance_in_x(head,-1,1);
                  head=advance_in_y(head,tail->y,0);
                  if (tail->y==head->y)
                     head=advance_in_z(head,tail->z,0);
                  else
                     head=advance_in_z(head,tail->z,1);
                } else
                    if (tail->x==columns-1) { // The tail is on the right
                       head=advance_in_x(head,columns-1,1);
                       head=advance_in_y(head,tail->y,0);
                       if (tail->y==head->y)
                         head=advance_in_z(head,tail->z,0);
                       else
                         head=advance_in_z(head,tail->z,1);
                    } else
                      if (tail->z==-1) { // The tail is in front 
                         head=advance_in_z(head,-1,1);
                         head=advance_in_y(head,(tail->y),0);
                         if (tail->y==head->y)
                            head=advance_in_x(head,tail->x,0);
                         else
                           head=advance_in_x(head,tail->x,1);
                      } else 
                         if (tail->z==depth-1) {
                          // The tail is at the back face
                           head=advance_in_z(head,depth-1,1);
                           head=advance_in_y(head,(tail->y),0);
                           if (tail->y==head->y)
                            head=advance_in_x(head,tail->x,0);
                           else
                            head=advance_in_x(head,tail->x,1);
                        } else {
                          // It is on the opposite side
                          minimum_option=0;
                          minimum_value=2*(depth-2)-(head->z)-(tail->z); //towards the end

                          this_value=(head->z)+(tail->z); // towards me
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=1;
                          }
             
                          this_value=(head->x)+(tail->x); // left
                          if (this_value<minimum_value){
                            minimum_value=this_value;
                            minimum_option=2;
                          }
             
                          this_value=2*(columns-2)-(head->x)-(tail->x); // right
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=3;
                          }
             
                          switch (minimum_option){
                              case 0: head=advance_in_z(head,depth-1,1);break;
                              case 1: head=advance_in_z(head,-1,1);break;
                              case 2: head=advance_in_x(head,-1,1);break;
                              case 3: head=advance_in_x(head,columns-1,1);break;
                         }
                         if (head->y==rows-1)
                           head=advance_in_y(head,-1,1);
                         else
                           head=advance_in_y(head,rows-1,1);
                         complete_loop(head,tail);
                      }
           } else
           if (head->z==-1 || head->z==depth-1) { // The head is at front or back faces
             if (tail->z==head->z) { // It is on the same side
                 head=advance_in_x(head,tail->x,0);
                 if (tail->x==head->x)
                     head=advance_in_y(head,tail->y,0);
                 else
                     head=advance_in_y(head,tail->y,1);

             } else
                if (tail->x==-1) { // The tail is on the left
                  head=advance_in_x(head,-1,1);
                  head=advance_in_z(head,tail->z,0);
                  if (tail->z==head->z)
                    head=advance_in_y(head,tail->y,0);
                  else
                    head=advance_in_y(head,tail->y,1);

                  head=advance_in_y(head,tail->y,1);
                } else
                    if (tail->x==columns-1) { // The tail is on the right
                       head=advance_in_x(head,columns-1,1);
                       head=advance_in_z(head,tail->z,0);
                       if (tail->z==head->z)
                         head=advance_in_y(head,tail->y,0);
                       else
                         head=advance_in_y(head,tail->y,1);

                    } else
                      if (tail->y==-1) { // The tail is on top
                         head=advance_in_y(head,-1,1);
                         head=advance_in_z(head,(tail->z),0);
                         if (tail->z==head->z)
                           head=advance_in_x(head,tail->x,0);
                         else
                           head=advance_in_x(head,tail->x,1);
                      } else 
                         if (tail->y==rows-1) {
                          // The tail is at the back face
                           head=advance_in_y(head,rows-1,1);
                           head=advance_in_z(head,(tail->z),0);
                           if (tail->z==head->z)
                             head=advance_in_x(head,tail->x,0);
                           else
                             head=advance_in_x(head,tail->x,1);
                        } else {
                          // It is on the opposite side
                          minimum_option=0;
                          minimum_value=2*(rows-2)-(head->y)-(tail->y); //towards the bottom

                          this_value=(head->y)+(tail->y); // towards the top
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=1;
                          }
             
                          this_value=(head->x)+(tail->x); // left
                          if (this_value<minimum_value){
                            minimum_value=this_value;
                            minimum_option=2;
                          }
             
                          this_value=2*(columns-2)-(head->x)-(tail->x); // right
                          if (this_value<minimum_value){
                             minimum_value=this_value;
                             minimum_option=3;
                          }
             
                          switch (minimum_option){
                              case 0: head=advance_in_y(head,rows-1,1);break;
                              case 1: head=advance_in_y(head,-1,1);break;
                              case 2: head=advance_in_x(head,-1,1);break;
                              case 3: head=advance_in_x(head,columns-1,1);break;
                          }
                          if (head->z==-1)
                            head=advance_in_z(head,-1,1);
                          else
                            head=advance_in_z(head,depth-1,1);
                      
                         complete_loop(head,tail);
                      }
           } 
}


/*************************************************************************************/
/*                           LOOP HANDLING - PROCESSING                              */  
/*************************************************************************************/

// Calculate the center of gravity of a set of nodes. In reality, it does not calculate
// the center of gravity. To obtain the center of gravity, it would be necessary to divide each
// coordinate by the number of nodes stored in the same structure (in field number_of_elements).

// It does not accept an empty list. If NULL ius provided, a runtime error will occur.

center_of_gravity_struct calculate_center_of_gravity(residue_node* aux) {
     
     center_of_gravity_struct center_of_gravity;
     residue_node* temp=aux;
     center_of_gravity.x=0;
     center_of_gravity.y=0;
     center_of_gravity.z=0;
     center_of_gravity.number_of_elements=0;
     while (aux->next!=temp){
           center_of_gravity.x+=aux->x;
           center_of_gravity.y+=aux->y;
           center_of_gravity.z+=aux->z;
           center_of_gravity.number_of_elements++;
           aux=aux->next;

     }
     center_of_gravity.x+=aux->x;
     center_of_gravity.y+=aux->y;
     center_of_gravity.z+=aux->z;
     center_of_gravity.number_of_elements++;
     return center_of_gravity;
}

// It returns if replacing the node original by replacement would reduce the area inside the path, following
// the criteria specified in the paper. Ignores the Z coordinate. This can be used with two nodes which have
// the same Z coordinate.
int reduces_areaXY(residue_node* replacement, residue_node* original, center_of_gravity_struct* center_of_gravity) {
    float temp;
    float area_original=fabs(original->x-((float)center_of_gravity->x)/(center_of_gravity->number_of_elements));
    area_original*=area_original;
    temp=fabs(original->y-((float)center_of_gravity->y)/(center_of_gravity->number_of_elements));
    area_original+=temp*temp;
    
    float area_replacement=fabs(replacement->x-((float)center_of_gravity->x)/(center_of_gravity->number_of_elements));
    area_replacement*= area_replacement;
    temp=fabs(replacement->y-((float)center_of_gravity->y)/(center_of_gravity->number_of_elements));
    area_replacement+=temp*temp;
 
    return (area_replacement < area_original);
    
}

// It returns if replacing the node original by replacement would reduce the area inside the path, following
// the criteria specified in the paper. Ignores the Y coordinate. This can be used with two nodes which have
// the same Y coordinate.
int reduces_areaXZ(residue_node* replacement, residue_node* original, center_of_gravity_struct* center_of_gravity) {
    float temp;
    float area_original=fabs(original->x-((float)center_of_gravity->x)/(center_of_gravity->number_of_elements));
    area_original*=area_original;
    temp=fabs(original->z-((float)center_of_gravity->z)/(center_of_gravity->number_of_elements));
    area_original+=temp*temp;
    
    float area_replacement=fabs(replacement->x-((float)center_of_gravity->x)/(center_of_gravity->number_of_elements));
    area_replacement*= area_replacement;
    temp=fabs(replacement->z-((float)center_of_gravity->z)/(center_of_gravity->number_of_elements));
    area_replacement+=temp*temp;
 
 
    return (area_replacement < area_original);
}

// It returns if replacing the node original by replacement would reduce the area inside the path, following
// the criteria specified in the paper. Ignores the X coordinate. This can be used with two nodes which have
// the same X coordinate.
int reduces_areaYZ(residue_node* replacement, residue_node* original, center_of_gravity_struct* center_of_gravity) {
    float temp;
    float area_original=fabs(original->y-((float)center_of_gravity->y)/(center_of_gravity->number_of_elements));
    area_original*=area_original;
    temp=fabs(original->z-((float)center_of_gravity->z)/(center_of_gravity->number_of_elements));
    area_original+=temp*temp;
    
    float area_replacement=fabs(replacement->y-((float)center_of_gravity->y)/(center_of_gravity->number_of_elements));
    area_replacement*= area_replacement;
    temp=fabs(replacement->z-((float)center_of_gravity->z)/(center_of_gravity->number_of_elements));
    area_replacement+=temp*temp;
 
    return (area_replacement < area_original);
}


// It solves all u-loops in a path, updating the center_of_gravity to be that of the
// remaining path
// process is a parameter by reference to indicate if any U loop has been solved
residue_node* process_u_loops(residue_node* residue, center_of_gravity_struct* center_of_gravity,int* processed) {
    int resta;
    int just_processed=1;
    residue_node *last_processed=residue;
    residue_node* current_residue=residue;
    residue_node* next_residue=residue->next->next->next;
    frontier* temp;
    *processed=0;
    // This operation is safe because the length of a loop must be >=4
    // To check if there is a U-shape we compare the coordinates of one pixel with the coordinates 
    // of the third pixel after it down the path. If the values in two axis are the same, then we
    // have a U-shape
        
    
    while ((current_residue->next!=next_residue) &&  (last_processed!=current_residue || just_processed)) { 
          
      // If the first condition is false, there are only two nodes in the list and thus there are no more U-shapes   
      // If the second condition is false, we have searched the entire list and there are no more U-shapes
      
        resta=(next_residue->value)-(current_residue->value);

        // The substraction of the position indicates the type of U-loop. This is called
        // resta
      
        if (resta==-1) {
                    temp=NULL;
                    // open below
                    if (current_residue->next->y < current_residue->y) {
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(current_residue->x)+columns*(current_residue->y)+twodsize*(current_residue->z);                                
                    }
                    else
                    // open above
                    if (current_residue->next->y > current_residue->y) {       
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(current_residue->x)+columns*(1+(current_residue->y))+twodsize*(current_residue->z);                                
                    } else 
                    //open back
                    if (current_residue->next->z < current_residue->z) {      
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(current_residue->x)+columns*(current_residue->y)+(columns*(rows-1))*(current_residue->z);
                    } else
                    //open front
                    if (current_residue->next->z > current_residue->z) {      
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(current_residue->x)+columns*(current_residue->y)+(columns*(rows-1))*(1+(current_residue->z));
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned)); // we change the flag of the appropriate node

                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;

                    //number_of_nodes-=2;
                    current_residue=current_residue->back->back; // The resolution of this U-loop may have created another
                    // new U-loop just before
                    last_processed=current_residue;
                    just_processed=1;
                    *processed=1;

           } else if (resta==1) {
                    temp=NULL;
                    // with respect to current_residue
                    if (current_residue->next->y < current_residue->y) {   // open below
                        if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                          temp=dfrontiers+(next_residue->x)+columns*(next_residue->y)+twodsize*(next_residue->z);                                
                    } else  
                    // open above
                    if (current_residue->next->y > current_residue->y) {      
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(next_residue->x)+columns*(1+(next_residue->y))+twodsize*(next_residue->z);                                
                    } else
                    //open back
                    if (current_residue->next->z < current_residue->z) {      
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(next_residue->x)+columns*(next_residue->y)+(columns*(rows-1))*(next_residue->z);
                    } else
                    //open front
                    if (current_residue->next->z > current_residue->z) {      
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(next_residue->x)+columns*(next_residue->y)+(columns*(rows-1))*(1+(next_residue->z));
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned));
                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;
                    
                    current_residue=current_residue->back->back;
                    last_processed=current_residue;
                    just_processed=1;
                    *processed=1;
            } else if (resta==-columns) { // next_residue is above current_residue
                    temp=NULL;

                    if (current_residue->next->x < current_residue->x) {  // open right
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(current_residue->x)+columns*(current_residue->y)+twodsize*(current_residue->z);
                    } else
                    if (current_residue->next->x > current_residue->x) {  // open left
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(current_residue->x)+1+columns*(current_residue->y)+twodsize*(current_residue->z);
                    } else  
                    if (current_residue->next->z < current_residue->z) { // open back
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(current_residue->x)+(columns-1)*(current_residue->y)+(columns-1)*rows*(current_residue->z);
                    } else  
                    if (current_residue->next->z > current_residue->z) {// open front
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(current_residue->x)+(columns-1)*(current_residue->y)+(columns-1)*rows*(1+(current_residue->z));
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned));
                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;
                    current_residue=current_residue->back->back;
                    last_processed=current_residue;
                    just_processed=1;                    
                    *processed=1;
            } else if (resta==columns) { // current_residue is above next_residue
                    temp=NULL;
                    // open left
                    if (current_residue->next->x < current_residue->x) {// open right
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(next_residue->x)+columns*(next_residue->y)+twodsize*(next_residue->z);
                    } else
                    if (current_residue->next->x > current_residue->x) {// open left
                      if ((current_residue->z!=(depth-1)) && (current_residue->z!=-1))
                        temp=dfrontiers+(next_residue->x)+1+columns*(next_residue->y)+twodsize*(next_residue->z);
                    } else  
                    if (current_residue->next->z < current_residue->z) {// open back
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(next_residue->x)+(columns-1)*(next_residue->y)+(columns-1)*rows*(next_residue->z);
                    } else  
                    if (current_residue->next->z > current_residue->z) {// open front
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(next_residue->x)+(columns-1)*(next_residue->y)+(columns-1)*rows*(1+(next_residue->z));
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned));
                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;
                    current_residue=current_residue->back->back;
                    last_processed=current_residue;
                    just_processed=1;                    
                    *processed=1;
             } else if (resta==-twodsize) { // next_residue is in front of current_residue
                    temp=NULL;
                    if (current_residue->next->x < current_residue->x) {// open right
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(current_residue->x)+columns*(current_residue->y)+columns*(rows-1)*(current_residue->z);
                    } else
                    if (current_residue->next->x > current_residue->x) {// open left
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(current_residue->x)+1+columns*(current_residue->y)+columns*(rows-1)*(current_residue->z);
                    } else                                    
                    if (current_residue->next->y < current_residue->y) {// open down
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(current_residue->x)+(columns-1)*(current_residue->y)+(columns-1)*rows*(current_residue->z);
                    } else  
                    if (current_residue->next->y > current_residue->y) {// open above
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(current_residue->x)+(columns-1)*(1+(current_residue->y))+(columns-1)*rows*(current_residue->z);
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned));
                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;
                    current_residue=current_residue->back->back;
                    last_processed=current_residue;
                    just_processed=1;                    
                    *processed=1;
             } else if (resta==twodsize) { // next_residue is behind current_residue
                    temp=NULL;
                    if (current_residue->next->x < current_residue->x) { // open righr
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(next_residue->x)+columns*(next_residue->y)+columns*(rows-1)*(next_residue->z);
                    } else
                    if (current_residue->next->x > current_residue->x) { // open left
                      if ((current_residue->y!=(rows-1)) && (current_residue->y!=-1))
                        temp=frontier_array+(next_residue->x)+1+columns*(next_residue->y)+columns*(rows-1)*(next_residue->z);
                    } else                                    
                    if (current_residue->next->y < current_residue->y) { // open down
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(next_residue->x)+(columns-1)*(next_residue->y)+(columns-1)*rows*(next_residue->z);
                    } else  
                    if (current_residue->next->y > current_residue->y) { // open above
                      if ((current_residue->x!=(columns-1)) && (current_residue->x!=-1))
                        temp=hfrontiers+(next_residue->x)+(columns-1)*(1+(next_residue->y))+(columns-1)*rows*(next_residue->z);
                    }
                    if (temp!=NULL)
                      temp->is_banned=(!(temp->is_banned));
                    (center_of_gravity->x)-=(next_residue->back->x)+(current_residue->next->x);
                    (center_of_gravity->y)-=(next_residue->back->y)+(current_residue->next->y);
                    (center_of_gravity->z)-=(next_residue->back->z)+(current_residue->next->z);
                    (center_of_gravity->number_of_elements)-=2;
                    delete(next_residue->back);
                    delete(current_residue->next);
                    next_residue->back=current_residue;
                    current_residue->next=next_residue;
                    current_residue=current_residue->back->back;
                    last_processed=current_residue;
                    just_processed=1;                    
                    *processed=1;
            } else {
                 // there is no U-loop at this node                                          
                 current_residue=current_residue->next;
                 next_residue=next_residue->next;
                 just_processed=0;
              }
                    
     } // end of while loop
     return last_processed;
} // end of function


// It solves a single l-loop in a path (the first found), updating the center_of_gravity to be that of the
// remaining path
residue_node* process_l_loop(residue_node* residue, center_of_gravity_struct* center_of_gravity) {
    int resta;
    int just_processed=1;
    residue_node *last_processed=residue;
    residue_node *possible_replacement=new residue_node;
    residue_node* current_residue=residue;
    residue_node* next_residue=residue->next->next;
    frontier* temp;
    // This operation is safe because the length of a loop must be >=4
    // To check if there is an L-shape we compare the coordinates of one pixel with the coordinates 
    // of the second pixel after it down the path. Depending on the location of these two pixels we
    // decide on whether there is an L-shape loop
        
    if (current_residue==next_residue) return NULL;
      // If this condition is true, there are only two nodes in the list and thus there are no more L-shapes 
    do{      
        resta=(next_residue->value)-(current_residue->value);
        // The substraction of the position indicates the type of U-loop. This is called
        // resta
      
        if (resta==columns+1) { // down right 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // |_  
                      possible_replacement->x=current_residue->x+1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)+=1-columns;
                         (center_of_gravity->x)++;
                         (center_of_gravity->y)--;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(1+(current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      
                      } 
                    } else { // ¨|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y+1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)-=1-columns;
                         (center_of_gravity->x)--;
                         (center_of_gravity->y)++;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(1+(current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        else
        if (resta==1-columns) { // up right
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // |¨
                      possible_replacement->x=current_residue->x+1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)+=columns+1;
                         (center_of_gravity->x)++;
                         (center_of_gravity->y)++;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(current_residue->y)+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                         
                      } 
                    } else { // _|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y-1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)-=columns+1;
                         (center_of_gravity->x)--;
                         (center_of_gravity->y)--;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(current_residue->y)+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      }                     } 
        }
        else
        if (resta==-columns-1) { // up left 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // ¨|
                      possible_replacement->x=current_residue->x-1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)+=columns-1;
                         (center_of_gravity->x)--;
                         (center_of_gravity->y)++;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // |_
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y-1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)-=columns-1;
                         (center_of_gravity->x)++;
                         (center_of_gravity->y)--;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
                    
        else
        if (resta==columns-1) { // down left 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // _|  
                      possible_replacement->x=current_residue->x-1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)+=-1-columns;
                         (center_of_gravity->x)--;
                         (center_of_gravity->y)--;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(1+(current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // ¨|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y+1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXY(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->y=possible_replacement->y;
                         (current_residue->next->value)-=-1-columns;
                         (center_of_gravity->x)++;
                         (center_of_gravity->y)++;
                         if ((current_residue->z!=-1) && (current_residue->z!=depth-1)){
                            temp=dfrontiers+twodsize*(current_residue->z)+columns*(1+(current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        // --------------------------------------------------------------------------------------
        
        else
        if (resta==twodsize+1) { // back right 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // /¨  
                      possible_replacement->x=current_residue->x+1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=1-twodsize;
                         (center_of_gravity->x)++;
                         (center_of_gravity->z)--;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(1+(current_residue->z))+columns*((current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         
                         return current_residue->back->back;
                      } 
                    } else { //   _/
                      possible_replacement->x=current_residue->x; 
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z+1;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=1-twodsize;
                         (center_of_gravity->x)--;
                         (center_of_gravity->z)++;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(1+(current_residue->z))+columns*((current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        else
        if (resta==1-twodsize) { // front right  
                    // To check the type of L 
                    if (current_residue->next->x == current_residue->x) { // /_
                      possible_replacement->x=current_residue->x+1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=1+twodsize;
                         (center_of_gravity->x)++;
                         (center_of_gravity->z)++;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(current_residue->z)+columns*((current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { //   ¨/
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z-1;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=1+twodsize;
                         (center_of_gravity->x)--;
                         (center_of_gravity->z)--;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(current_residue->z)+columns*((current_residue->y))+(current_residue->x)+1;
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        
        else
        if (resta==-twodsize-1) { // front left 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // _/
                      possible_replacement->x=current_residue->x-1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=twodsize-1;
                         (center_of_gravity->x)--;
                         (center_of_gravity->z)++;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(current_residue->z)+columns*((current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // /¨
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z-1;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=twodsize-1;
                         (center_of_gravity->x)++;
                         (center_of_gravity->z)--;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(current_residue->z)+columns*((current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
                    
        else
        if (resta==twodsize-1) { // back left 
                    // To check the type of L
                    if (current_residue->next->x == current_residue->x) { // ¨/
                      possible_replacement->x=current_residue->x-1;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=1+twodsize;
                         (center_of_gravity->x)--;
                         (center_of_gravity->z)--;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(1+(current_residue->z))+columns*((current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         
                         return current_residue->back->back;
                      } 
                    } else { // /_
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z+1;
                      if (reduces_areaXZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->x=possible_replacement->x;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=1+twodsize;
                         (center_of_gravity->x)++;
                         (center_of_gravity->z)++;
                         if ((current_residue->y!=-1) && (current_residue->y!=rows-1)){
                            temp=frontier_array+((rows-1)*columns)*(1+(current_residue->z))+columns*((current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        
        // --------------------------------------------------------------------------------------
        
        else
        if (resta==twodsize+columns) { // back down 
                    // To check the type of L
                    if (current_residue->next->y == current_residue->y) { // /|  
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y+1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=columns-twodsize;
                         (center_of_gravity->y)++;
                         (center_of_gravity->z)--;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(1+(current_residue->z))+(columns-1)*((current_residue->y)+1)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // |/
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z+1;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=columns-twodsize;
                         (center_of_gravity->y)--;
                         (center_of_gravity->z)++;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(1+(current_residue->z))+(columns-1)*((current_residue->y)+1)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        else
        if (resta==twodsize-columns) { // back up
                    // To check the type of L
                    if (current_residue->next->y == current_residue->y) { // .|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y-1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=-twodsize-columns;
                         (center_of_gravity->y)--;
                         (center_of_gravity->z)--;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(1+(current_residue->z))+(columns-1)*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // |´
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z+1;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=-twodsize-columns;
                         (center_of_gravity->y)++;
                         (center_of_gravity->z)++;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(1+(current_residue->z))+(columns-1)*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        
        else
        if (resta==-twodsize-columns) { // front up
                    // To check the type of L
                    if (current_residue->next->y == current_residue->y) { // |.
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y-1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=twodsize-columns;
                         (center_of_gravity->y)--;
                         (center_of_gravity->z)++;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(current_residue->z)+(columns-1)*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } else { // .|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z-1;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=twodsize-columns;
                         (center_of_gravity->y)++;
                         (center_of_gravity->z)--;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(current_residue->z)+(columns-1)*(current_residue->y)+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        else                   
        if (resta==-twodsize+columns) { // front down
                    // To check the type of L
                    if (current_residue->next->y == current_residue->y) { // |´  
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y+1;
                      possible_replacement->z=current_residue->z;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)+=twodsize+columns;
                         (center_of_gravity->y)++;
                         (center_of_gravity->z)++;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(current_residue->z)+(columns-1)*(1+(current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
        
                    } else { // .|
                      possible_replacement->x=current_residue->x;
                      possible_replacement->y=current_residue->y;
                      possible_replacement->z=current_residue->z-1;
                      if (reduces_areaYZ(possible_replacement, current_residue->next, center_of_gravity)) {
                         current_residue->next->y=possible_replacement->y;
                         current_residue->next->z=possible_replacement->z;
                         (current_residue->next->value)-=twodsize+columns;
                         (center_of_gravity->y)--;
                         (center_of_gravity->z)--;
                         if ((current_residue->x!=-1) && (current_residue->x!=columns-1)){
                            temp=hfrontiers+((columns-1)*rows)*(current_residue->z)+(columns-1)*(1+(current_residue->y))+(current_residue->x);
                            temp->is_banned=!(temp->is_banned);
                         }                           
                         return current_residue->back->back;
                      } 
                    } 
        }
        
        current_residue=current_residue->next;
        next_residue=next_residue->next;
     } while (current_residue!=residue); // end of while loop

     // All the L loops found are no good (do not reduce the surface). If there are no U loops and this happens,
     // it is because all L loops leave the total surface unaltered. We have to make a decision. We alter the
     // center of gravity slightly
     // This could have been a random increment, but to make it the algorithm deterministic we always add a 
     // small increment to each axis.
    
     (center_of_gravity->x)=(center_of_gravity->x)+0.001;
     (center_of_gravity->y)=(center_of_gravity->y)+0.001;
     (center_of_gravity->z)=(center_of_gravity->z)+0.001;
     return process_l_loop(residue, center_of_gravity);
} // end of function



// It processes a loop, according to Hussein's thesis. This operation consists of alternating between solving U-loops
// and L-Loops, until the path is reduced to two nodes.
void process_loop(residue_node* residue) {

   int u_loop_solved;
   center_of_gravity_struct c;
   if (residue!=NULL)
     c=calculate_center_of_gravity(residue);
   while (residue!=NULL) {
         u_loop_solved=0;
         residue=process_u_loops(residue,&c,&u_loop_solved);
         if (!u_loop_solved)      {
           // If no U loop has been solved, we move the pointer forward so that the last L-Loop solved
           // is the last one to be checked
           residue=residue->next->next->next;
         }
         residue=process_l_loop(residue,&c);
   } // end of while
}


// Builds the path and processes it
void build_and_process_loop(residue_values* residue) {
         residue_node* forward;     
         residue_node* backward;
         residue_node* temp;
         residue_node* new_node;

         forward=build_loop_forward(residue); // build a path following output lines
         if (forward==NULL) // this box doe not contain any residues          
           return;
         backward=build_loop_backward(residue); // build a path following input lines

         if (backward!=NULL) { // loop breaks at the border (otherwise 
            // build_loop_forward would have taken the node)
            temp=backward; // to store the last element of the list
            while (temp->next!=NULL)
              temp=temp->next;

            temp->next=forward->next;  // forward cannot be NULL. Otherwise
            // we would not have entered the if statement
            forward->next->back=temp;

            // first input and output correspond to the same cube
            // we delete one of them at the same time as the lists are
            // joined
            delete(forward);

             while (temp->next!=NULL)
                temp=temp->next;  // Move to the last element of the list
                               // we start at temp to make it faster
             // At the end of the loop, tail contains one border and
             // temp the other
             complete_loop(temp,backward); // It completes the loop if necessary
         } else {
             temp=forward;
             backward=forward;
         } 
         
         // Take temp to the end of the list
         while (temp->next!=NULL)
           temp=temp->next; 
            
         // In whichever case, temp points to the tail and backward to the first element
         // Make it a round list
         if (backward!=NULL) {
             backward->back=temp;
             temp->next=backward;
         }                
         
         // Now that the loop has been built. We can process it.
         process_loop(backward); 
         
}   


     


// The loop that calls the function that build and process all loops
void build_and_process_all_loops() {
   residue_node* last=NULL;
   residue_node* temp=NULL;
   residue_values* current_residue=residues;
   // Try all possible residues as seeds and build paths
   for (int d=0;d<depth-1;d++) {
     for (int r=0;r<rows-1;r++) {
       for (int c=0;c<columns-1;c++) {
           build_and_process_loop(current_residue);
          current_residue++;
       }
       current_residue++;
     }
     current_residue+=columns;
   }
}              


/*************************************************************************************/
/*                         INCREMENT HANDLING FUNCTIONS                              */  
/*************************************************************************************/
group* getAndSortGroup(group* the_group) {
   if (the_group->previous_group==NULL)
	   return the_group;
   group* old_previous=the_group->previous_group;
   group* result=getAndSortGroup(the_group->previous_group);
   the_group->previous_group=result;
   the_group->current_increment+=old_previous->current_increment;
   return result;
}



// Resolves all values in a group with respect to its first element. At the end
// of the process, each pixel will hold in the field current_increment the increment
// (the number of 2pi jumps which should be added) with respect to the first element 
// of the group (the head). This means that the current value of the pixel would be
// its value plus 2pi multiplied by the addition of the values stored in the field 
// current_increment of the pixel plus the value stored in the field current_increment
// of the pixel pointed by the pointer in previous_pixel. Note that at the end of
// the execution the field previous_pixel in all pixels in the group will point
// to the head.
void solve_image(group* gr)
{
   group* first_element=gr; // gr->current_increment debe ser siempre 
   gr=gr->next_pixel;
   while (gr!=NULL)
   {
     gr->current_increment+=gr->previous_group->current_increment;
     gr->previous_group=first_element;
	 gr=gr->next_pixel;
   }
//  printf("number of groups processed: %d\n",i);

}

void calculate_final_pixel_values()
{
  group *p=array;
  unsigned char* m=mask;
  if (m!=NULL) {
      while (m==0) {
            m++;
            p++;
      }      
      if (m-mask==number_of_pixels){
       printf("\n ***** No output produced: the mask is all zeros ****\n");
       exit(0);
      } 
  }
   while (p->previous_group!=NULL) {
       p=p->previous_group;
   }
  
  // Once this loop finishes we have the head of the group
  solve_image(p);
}


// Calculates the real value of a pixel accumulating the increments until the 
// root is reached. 
// This value is obtained // combining the original pixel value plus a number of 2pi 
// jumps. The number of 2pi // jumps to add to the pixel is given by its current_increment 
// plus the increment of the pixel pointed by the field previous_pixel. Note that this
// increment also needs to be calculated in the same manner (it is a recursive definition)
float calculate_pixel_value(group *gr)
{
	int increment=0;
	group* temp=gr;
	while (temp!=NULL)
	{
        increment+=temp->current_increment;
		temp=temp->previous_group;
	}

	return (gr->value)+(twopi*increment);
}

/*************************************************************************************/
/*                   FUNCTIONS ABOUT CALCULATING PIXEL DISTANCES                     */  
/*************************************************************************************/

// Given two values, it returns the potential minimum distance between them after
// a local unwrap operation (in dist). It also provides the number of 2pi steps which would
// be required to add to the second element so that this distance is obtained (in steps)
void calculate_complete_distance(float number1,float number2, float* dist, int* steps)
{

	float diff=(number1-number2);
    float difference=(diff)/twopi;
    *steps=(int)(floor)(difference+0.5);
    *dist=diff-(*steps)*twopi;
}


// The same as above but it does not return the number of 2pi steps
float calculate_distance(float number1,float number2)
{
    int steps;
	float diff=(number1-number2);
    float difference=(diff)/twopi;
    steps=(int)(floor)(difference+0.5);
    return diff-(steps)*twopi;
}


// It only calculates the number of 2pi steps
int calculate_steps(float number1,float number2)
{
    float difference=(number1-number2)/twopi;
    return (int)floor(difference+0.5);
}





/*************************************************************************************/
/*                    FUNCTIONS TO GENERATE THE STRUCTURE                            */  
/*************************************************************************************/
void calculate_derivatives()
{
   
   // We place current_pixel in the first row first column of the image
   group* current_pixel=array;
   residue_values* current_residue=residues;
   group* neighbour_down=array+columns;
   group* neighbour_in_depth=array+twodsize;
   distances* current_derivative=derivatives;
   float current_pixel_value;
   float neighbour_in_depth_pixel_value;
   float neighbour_down_pixel_value;
   float temp;
   float dummy;
   if (debug)
     printf("Starting computation of differences...");
      
   current_pixel_value=current_pixel->value;
   neighbour_down_pixel_value=neighbour_down->value; // The value of the pixel below
   neighbour_in_depth_pixel_value=neighbour_in_depth->value; // The value of the next pixel in 3rd dimension
   int steps;      
   for (int d=0;d<depth-1;d++) {
     for (int r=0;r<rows-1;r++) {
       for (int c=0;c<columns-1;c++) {
       calculate_complete_distance(current_pixel_value,neighbour_down_pixel_value,
              &(current_derivative->vertical_distance), &steps);

          current_derivative->vertical_steps=steps;
          current_residue->residue_front_sum=-steps;
          current_residue->residue_side_sum=steps;
          if (c!=0) 
            (current_residue-1)->residue_front_sum+=steps;
          if (d!=0)
            (current_residue-twodsize)->residue_side_sum-=steps;
          
          calculate_complete_distance(current_pixel_value,neighbour_in_depth_pixel_value,
              &(current_derivative->in_depth_distance), &steps);

          current_derivative->in_depth_steps=steps;
          current_residue->residue_top_sum=steps;
          current_residue->residue_side_sum-=steps;
          
          if (c!=0)
            (current_residue-1)->residue_top_sum-=steps;
          if (r!=0)
            (current_residue-columns)->residue_side_sum+=steps;  
             
          
          temp=(neighbour_in_depth+columns)->value;
          current_derivative->diagonal_right_side_distance=calculate_distance(current_pixel_value,temp);
              
          current_derivative->diagonal_left_side_distance=calculate_distance(neighbour_in_depth_pixel_value, neighbour_down_pixel_value);
      
          neighbour_in_depth++;
      
          current_derivative->diagonal_right_above_distance=calculate_distance(current_pixel_value, neighbour_in_depth->value);
      
          current_pixel++;
      
          current_derivative->diagonal_left_distance=calculate_distance(current_pixel->value, neighbour_down_pixel_value);
      
          calculate_complete_distance(current_pixel_value,current_pixel->value,
              &(current_derivative->horizontal_distance), &steps);

          current_derivative->horizontal_steps=steps;
          current_residue->residue_front_sum+=steps;
          current_residue->residue_top_sum-=steps;
          if (r!=0)
            (current_residue-columns)->residue_front_sum-=steps;
          
          if (d!=0)
            (current_residue-twodsize)->residue_top_sum+=steps;
                    
        
          current_derivative->opposite_corner_from_bottom=calculate_distance(neighbour_down_pixel_value, neighbour_in_depth->value);
          
          neighbour_down++;
      
          neighbour_down_pixel_value=neighbour_down->value; // The value of the pixel below
      
          current_derivative->opposite_corner_from_opposite=calculate_distance(neighbour_down_pixel_value, neighbour_in_depth_pixel_value);
 
          current_derivative->diagonal_right_distance=calculate_distance(current_pixel_value, neighbour_down_pixel_value);
 
          current_derivative->opposite_corner_from_me=calculate_distance(current_pixel_value, (neighbour_in_depth+columns)->value);
 
          current_pixel_value=current_pixel->value;

          current_derivative->opposite_corner_from_right=calculate_distance(current_pixel_value, temp);
          
          current_derivative->diagonal_left_above_distance=calculate_distance(current_pixel_value, neighbour_in_depth_pixel_value);
      
          neighbour_in_depth_pixel_value=neighbour_in_depth->value; // The value of the pixel below
          
          current_derivative++;
          
          current_residue++;
        }

        steps=calculate_steps(current_pixel_value,neighbour_down_pixel_value);

        current_residue->residue_side_sum=steps;
        (current_residue-1)->residue_front_sum+=steps;
        if (d!=0)
         (current_residue-twodsize)->residue_side_sum-=steps;
          
        steps=calculate_steps(current_pixel_value,neighbour_in_depth_pixel_value);
        
        current_residue->residue_side_sum-=steps;
          
        (current_residue-1)->residue_top_sum-=steps;
        if (r!=0)
          (current_residue-columns)->residue_side_sum+=steps;  
             
        neighbour_down++;   
        current_pixel++; 
          
        current_residue++;
        neighbour_in_depth++;
        current_pixel_value=current_pixel->value;
        neighbour_down_pixel_value=neighbour_down->value; // The value of the pixel below
        neighbour_in_depth_pixel_value=neighbour_in_depth->value; // The value of the pixel below
      }
      // For the last row of each plane:
      for (int c=0;c<columns-1;c++) {
          steps=calculate_steps(current_pixel->value,(current_pixel+twodsize)->value);
          current_residue->residue_top_sum=steps;
          
          if (c!=0)
            (current_residue-1)->residue_top_sum-=steps;
          (current_residue-columns)->residue_side_sum+=steps;  
           
          steps=calculate_steps(current_pixel->value,(current_pixel+1)->value);
          current_residue->residue_top_sum-=steps;
          (current_residue-columns)->residue_front_sum-=steps;
          if (d!=0)
            (current_residue-twodsize)->residue_top_sum+=steps;
                    
          current_pixel++;  
          current_residue++;

      }
      steps=calculate_steps(current_pixel->value,(current_pixel+twodsize)->value);
      (current_residue-1)->residue_top_sum-=steps;
      (current_residue-columns)->residue_side_sum+=steps;  

      neighbour_down+=columns;   
      current_pixel++;   
      current_residue++;
      neighbour_in_depth+=columns;
      current_pixel_value=current_pixel->value;
      neighbour_down_pixel_value=neighbour_down->value; // The value of the pixel below
      neighbour_in_depth_pixel_value=neighbour_in_depth->value; // The value of the pixel below
      
   }
   
   // Process the back layer
   
   for (int r=0;r<rows-1;r++) {
     for (int c=0;c<columns-1;c++) {

          steps=calculate_steps(current_pixel->value,(current_pixel+columns)->value);

          current_residue->residue_front_sum=-steps;
          if (c!=0) 
            (current_residue-1)->residue_front_sum+=steps;
          (current_residue-twodsize)->residue_side_sum-=steps;
          
          steps=calculate_steps(current_pixel->value,(current_pixel+1)->value);
          current_residue->residue_front_sum+=steps;
          if (r!=0)
            (current_residue-columns)->residue_front_sum-=steps;
          
          (current_residue-twodsize)->residue_top_sum+=steps;
          current_pixel++;
          current_residue++;
     }

     steps=calculate_steps(current_pixel->value,(current_pixel+columns)->value);

     (current_residue-1)->residue_front_sum+=steps;
     (current_residue-twodsize)->residue_side_sum-=steps;
             
     current_pixel++;   
     current_residue++;
  }
     
  // For the last row of the last plane
     for (int c=0;c<columns-1;c++) {

          steps=calculate_steps(current_pixel->value,(current_pixel+1)->value);
          (current_residue-columns)->residue_front_sum-=steps;
          
          (current_residue-twodsize)->residue_top_sum+=steps;
          current_pixel++;
          current_residue++;
     }

   if (debug)
     printf("  finished\n");

}
// Calculates the reliability of each pixel, according to the second derivatives.
// The second derivatives in all directions are calculated and the reliability is
// calculated as the inverse of the addition of the squares

void calculate_pixel_reliability()
{
   
   int i,j;
   if (debug)
      printf("Starting computation of reliability values...\n");
   calculate_derivatives();
   // We place current_pixel in the first row first column of the image
   group* current_pixel=array;
   //group* neighbour_down=array+columns;
   /* calculate the distances */
      distances* current_derivative=derivatives;
      float current_pixel_value;
      float neighbour_down_pixel_value;
      current_pixel_value=current_pixel->value;
      
   
    /* The reliability values are calculated */
      // In the distances we have not calculated the last pixels of each row
   
      // Now we have to store the realibility for each pixel
      // Inverse values are stored: higher values imply a lower reliability 
      // Pixels in the borders of the image are assigned very high values (very low reliability)
   
      // The reliability values assigned are the addition of the absolute values of the second derivatives in 
      // each direction
   
      
      if (debug)
         printf("  Top layer...\n");
   
      /* top layer */
      for (int d=0;d<depth;d++) {
        current_pixel=array+(d*twodsize);        
        for (int c=0;c<columns;c++) {
           current_pixel->reliability=MIN_RELIABILITY;
           current_pixel++;
         }
      }

      if (debug)
         printf("  Bottom layer...\n");
      
      /* bottom layer */
      for (int d=1;d<=depth;d++) {
        current_pixel=array+(d*twodsize)-columns;        
        for (int c=0;c<columns;c++) {
           current_pixel->reliability=MIN_RELIABILITY;
           current_pixel++;
         }
      }

      if (debug)
         printf("  Front layer...\n");
      
      /* front layer */
      current_pixel=array+columns;
      for (int r=1;r<rows-1;r++) { // La de arriba y la de abajo ya están
        for (int c=0;c<columns;c++) {
           current_pixel->reliability=MIN_RELIABILITY;
           current_pixel++;
         }
      }
      
      if (debug)
         printf("  Back layer...\n");

      /* back layer */
      current_pixel=array+number_of_pixels-twodsize+columns;
      for (int r=1;r<rows-1;r++) { // La de arriba y la de abajo ya están
        for (int c=0;c<columns;c++) {
           current_pixel->reliability=MIN_RELIABILITY;
           current_pixel++;
         }
      }

      // Right and left layers are processed with the rest of the image
     if (debug)
         printf("  Rest of the image...\n");

      /* rest of the image */
         float accumulated_derivative;
         float partial_derivative;
         current_derivative;
         distances* temp_derivative;
         // Rest of the image except last row
      int temp=(columns-1)*(rows-1);

      for (int d=1;d<depth-1;d++){
         current_pixel=array+(d*twodsize)+columns; // Place current_pixel                                                    
         
         current_derivative=derivatives+(d*temp)+columns; // The one that
         

         // corresponds to current_pixel+1
         
         for (int r=1;r<rows-1;r++) { // The first and last rows are excluded
           // First and last pixel of each column have infinite values
           /* First column of the row */
           current_pixel->reliability=MIN_RELIABILITY; // El primero de la fila
           current_pixel++;
           
           for (int c=1;c<columns-1;c++) {
               // FRONT CUBE
               /* right diagonal */
               accumulated_derivative=current_derivative->diagonal_right_distance-
                 (current_derivative-columns)->diagonal_right_distance;
               accumulated_derivative*=accumulated_derivative;
               
               /* vertical */
               partial_derivative=current_derivative->vertical_distance-
                   (current_derivative-columns+1)->vertical_distance;
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
                 
               /* left diagonal */
               partial_derivative=(current_derivative-columns+1)->diagonal_left_distance-
                 (current_derivative-1)->diagonal_left_distance;
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
                 
               /* horizontal */
               partial_derivative=(current_derivative)->horizontal_distance-
                 (current_derivative-1)->horizontal_distance;
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;


               // SIDE CUBE
               /* right diagonal */
               partial_derivative=current_derivative->diagonal_right_side_distance-
                   (current_derivative-temp-columns+1)->diagonal_right_side_distance;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
               /* vertical */
               partial_derivative=current_derivative->in_depth_distance-
                   (current_derivative-temp)->in_depth_distance;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
                 
               /* left diagonal */
               partial_derivative=(current_derivative-columns+1)->diagonal_left_side_distance-
                 (current_derivative-temp)->diagonal_left_side_distance;
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
                 

               
               // TOP CUBE
               /* right diagonal */
               partial_derivative=current_derivative->diagonal_right_above_distance-
                   (current_derivative-temp-1)->diagonal_right_above_distance;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
               /* left diagonal */
               partial_derivative=(current_derivative-1)->diagonal_left_above_distance-
                 (current_derivative-temp)->diagonal_left_above_distance;
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;

               // LAST DIAGONALS
               //opposite_corner_from_me;
               partial_derivative=current_derivative->opposite_corner_from_me-
                   (current_derivative-temp-columns)->opposite_corner_from_me;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
               //opposite_corner_from_right;
               partial_derivative=(current_derivative-temp-columns+1)->opposite_corner_from_right-
                   (current_derivative-1)->opposite_corner_from_right;
               
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
               //opposite_corner_from_bottom;
               partial_derivative=(current_derivative-columns+1)->opposite_corner_from_bottom-
                   (current_derivative-temp-1)->opposite_corner_from_bottom;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
               // opposite_corner_from_opposite;
               partial_derivative=(current_derivative-columns)->opposite_corner_from_opposite-
                   (current_derivative-temp)->opposite_corner_from_opposite;
                   
               partial_derivative*=partial_derivative;
               accumulated_derivative+=partial_derivative;
               
                 /* assignement */
                 current_pixel->reliability=accumulated_derivative;
                 current_pixel++;
                 current_derivative++;
           }
           
           /* Last column of the row */
           current_pixel->reliability=MIN_RELIABILITY;
           current_pixel++;
          current_derivative++;
         }
      
    }      

   if (debug)
     printf("finished.\n");
      
   
}


// Initial creation of borders. Initially, it exists a border between any two adjacent 
// pixels in the image. Each border has a reliability which depends upon the reliability
// of the pixels at both sides.

void build_frontiers()
{
   
   int i,j;
   if (debug) 
     printf("Building frontiers...\n");
   frontier* new_frontier=frontier_array;
   group* temp_group_current;
   group* temp_group_compareTo;
   
   
   float distance;
   
   /* First we create the borders of each pixels with the one inmediately below */
      distances* temp_derivatives=derivatives;
      
      temp_group_current=array;  // the current pixel
      temp_group_compareTo=array+columns;  // the pixel immediately below
   
      if (debug) 
       printf("  Vertical frontiers...\n");

      for (int d=0;d<depth-1;d++) {
   
        for (int r=0;r<rows-1;r++) {
            for (int c=0;c<columns-1;c++) { // the last column is processed after the loop
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
 
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=temp_derivatives->vertical_steps;
              new_frontier++;
			
              temp_group_current++;
              temp_group_compareTo++;
              temp_derivatives++; // *** posible causa de error. revisar en caso
              // de que aparezcan falsos wraps en el resultado
            }
            /* We create a new border, with the appropriate pixels */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
         
            /* The distance is given by the addition of the squares of the second derivatives in all directions */
            // The vertical steps of the last column were not calculated. They need to be calculated now
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
         
            new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
         
            new_frontier++;
         
            temp_group_current++;
            temp_group_compareTo++;
         } // end processing a row
         // lets advance to the next 2D plane:
         // Last row does not need to be process (there are no borders to add)
         temp_group_current+=columns; // one entire row
         temp_group_compareTo+=columns;
         // temp_derivatives must be on place
                 
      }
      //The back layer
        for (int r=0;r<rows-1;r++) {
            for (int c=0;c<columns;c++) { // the last column is processed after the loop
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
            
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
            
              new_frontier++;
              
              temp_group_current++;
              temp_group_compareTo++;
              
            }
        }

   if (debug)
     printf("  Horizontal frontiers...\n");

      /* The borders of each pixels with the one inmediately to the right */
      temp_derivatives=derivatives;
      temp_group_current=array; // the current pixel
      temp_group_compareTo=array+1; // the pixel to the right
      for (int d=0;d<depth-1;d++) {
        for (int i=0;i<rows-1;i++) {
          for (int j=0;j<columns-1;j++) {
            /* Update the values of the groups at both side of the border. Since there are only two pixels separated
            by a border, each group will be composed of a single pixel */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
          
            /* Calculate the distance */
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
            
            new_frontier->steps=temp_derivatives->horizontal_steps;
          
            /* Insert the border in the heap */
            new_frontier++;
            temp_group_current++;
            temp_group_compareTo++;
            temp_derivatives++;
         } // a row has been processed
         temp_group_current++;
         temp_group_compareTo++;
       } // a 2D plane has been processed except for the last row
       // There are no pre-calculated distances for the last row
       for (j=0;j<columns-1;j++) {
          /* We update the pointers with the groups, which are both composed of a single pixel */
          new_frontier->group1=temp_group_current;
          new_frontier->group2=temp_group_compareTo;
        
          /* Update the distance between the pixels at both side of the border */
          distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          //distance*=distance;
          new_frontier->distance=distance;
          new_frontier->distance_backup=0;
          new_frontier->is_banned=0;
          
          new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
        
          /* Insert the new border in the heap*/
          new_frontier++;
          temp_group_current++;
          temp_group_compareTo++;
          //temp_derivatives++;
       } // Now the entire 2D plane has been processed
       // Update pointers
       temp_group_current++;
       temp_group_compareTo++;
   } // All 2D planes have been processed
   
     //The back layer
        for (int r=0;r<rows;r++) {
            for (int c=0;c<columns-1;c++) { 
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
            
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
            
              new_frontier++;
			
              temp_group_current++;
              temp_group_compareTo++;
              
            }
            temp_group_current++;
            temp_group_compareTo++;

        }
   
   
   
   if (debug)
     printf("  Z frontiers...\n");
   
      /* The borders of each pixels with the one in the inmediate 2D plane */
      temp_derivatives=derivatives;
      temp_group_current=array; // the current pixel
      temp_group_compareTo=array+twodsize; // the pixel to the right
      for (int d=0;d<depth-1;d++) {
        for (int i=0;i<rows-1;i++) {
          for (int j=0;j<columns-1;j++) {
            /* Update the values of the groups at both side of the border. Since there are only two pixels separated
            by a border, each group will be composed of a single pixel */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
          
            /* Calculate the distance */
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
            
            new_frontier->steps=temp_derivatives->in_depth_steps;
          
            /* Insert the border in the heap */
            new_frontier++;
            temp_group_current++;
            temp_group_compareTo++;
            temp_derivatives++;
         } // a row has been processed except for the last column
         /* We create a new border, with the appropriate pixels */
         new_frontier->group1=temp_group_current;
         new_frontier->group2=temp_group_compareTo;
         
         /* The distance is given by the addition of the squares of the second derivatives in all directions */
         // The vertical steps of the last column were not calculated. They need to be calculated now
         distance=temp_group_current->reliability+temp_group_compareTo->reliability;
         //distance*=distance;
         new_frontier->distance=distance;
         new_frontier->distance_backup=0;
         new_frontier->is_banned=0;
         
         new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
         
         /* insert the border in the heap */
         new_frontier++;
        
         temp_group_current++;
         temp_group_compareTo++;
         
         //temp_group_current++;
         //temp_group_compareTo++;
       } // a 2D plane has been processed except for the last row
       // There are no pre-calculated distances for the last row
       for (j=0;j<columns;j++) {
          /* We update the pointers with the groups, which are both composed of a single pixel */
          new_frontier->group1=temp_group_current;
          new_frontier->group2=temp_group_compareTo;
        
          /* Update the distance between the pixels at both side of the border */
          distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          //distance*=distance;
          new_frontier->distance=distance;
          new_frontier->distance_backup=0;
          new_frontier->is_banned=0;
          
          new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
        
          /* Insert the new border in the heap*/
          new_frontier++;
          temp_group_current++;
          temp_group_compareTo++;
          //temp_derivatives++;
       } // Now the entire 2D plane has been processed
       // Update pointers
       //temp_group_current++;
       //temp_group_compareTo++;
   } // All 2D planes have been processed
   
   if (debug) {
     printf("A total of %d frontiers have been created\n",new_frontier-frontier_array);
   }

}

// Initial creation of borders. Initially, it exists a border between any two adjacent 
// pixels in the image. Each border has a reliability which depends upon the reliability
// of the pixels at both sides.

void build_frontiers(unsigned char* mask,float* qualities)
{
   
   int i,j;
   if (debug) 
     printf("Building frontiers...\n");
   frontier* new_frontier=frontier_array;
   group* temp_group_current;
   group* temp_group_compareTo;
   
   float distance;
   
   /* First we create the borders of each pixels with the one inmediately below */
      distances* temp_derivatives=derivatives;
      
      temp_group_current=array;  // the current pixel
      temp_group_compareTo=array+columns;  // the pixel immediately below
      
      
      if (debug) 
       printf("  Vertical frontiers...\n");

      for (int d=0;d<depth-1;d++) {
   
        for (int r=0;r<rows-1;r++) {
            for (int c=0;c<columns-1;c++) { // the last column is processed after the loop
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
 
              if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
                distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              else
                 distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
              if (mask!=NULL) {
                 if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                      distance=FLT_MAX;
              }
                   
                
                   
              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=temp_derivatives->vertical_steps;
              new_frontier++;
			
              temp_group_current++;
              temp_group_compareTo++;
              temp_derivatives++; // *** posible causa de error. revisar en caso
              // de que aparezcan falsos wraps en el resultado
            }
            /* We create a new border, with the appropriate pixels */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
         
            /* The distance is given by the addition of the squares of the second derivatives in all directions */
            // The vertical steps of the last column were not calculated. They need to be calculated now
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            else
              distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
            if (mask!=NULL) {
               if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                    distance=FLT_MAX;
            }

            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
         
            new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
         
            new_frontier++;
         
            temp_group_current++;
            temp_group_compareTo++;
         } // end processing a row
         // lets advance to the next 2D plane:
         // Last row does not need to be process (there are no borders to add)
         temp_group_current+=columns; // one entire row
         temp_group_compareTo+=columns;
         // temp_derivatives must be on place
                 
      }
      //The back layer
        for (int r=0;r<rows-1;r++) {
            for (int c=0;c<columns;c++) { // the last column is processed after the loop
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
            
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;

              if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
                distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              else
                 distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
              if (mask!=NULL) {
                 if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                      distance=FLT_MAX;
              }

              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
            
              new_frontier++;
              
              temp_group_current++;
              temp_group_compareTo++;
              
            }
        }

   if (debug)
     printf("  Horizontal frontiers...\n");

      /* The borders of each pixels with the one inmediately to the right */
      temp_derivatives=derivatives;
      temp_group_current=array; // the current pixel
      temp_group_compareTo=array+1; // the pixel to the right
      for (int d=0;d<depth-1;d++) {
        for (int i=0;i<rows-1;i++) {
          for (int j=0;j<columns-1;j++) {
            /* Update the values of the groups at both side of the border. Since there are only two pixels separated
            by a border, each group will be composed of a single pixel */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
          
            /* Calculate the distance */
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;

            if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            else
               distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
            if (mask!=NULL) {
               if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                    distance=FLT_MAX;
            }

            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
            
            new_frontier->steps=temp_derivatives->horizontal_steps;
          
            /* Insert the border in the heap */
            new_frontier++;
            temp_group_current++;
            temp_group_compareTo++;
            temp_derivatives++;
         } // a row has been processed
         temp_group_current++;
         temp_group_compareTo++;
       } // a 2D plane has been processed except for the last row
       // There are no pre-calculated distances for the last row
       for (j=0;j<columns-1;j++) {
          /* We update the pointers with the groups, which are both composed of a single pixel */
          new_frontier->group1=temp_group_current;
          new_frontier->group2=temp_group_compareTo;
        
          /* Update the distance between the pixels at both side of the border */
          distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          if (qualities==NULL)
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          else
            distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
          if (mask!=NULL) {
             if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                 distance=FLT_MAX;
          }

          //distance*=distance;
          new_frontier->distance=distance;
          new_frontier->distance_backup=0;
          new_frontier->is_banned=0;
          
          new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
        
          /* Insert the new border in the heap*/
          new_frontier++;
          temp_group_current++;
          temp_group_compareTo++;
          //temp_derivatives++;
       } // Now the entire 2D plane has been processed
       // Update pointers
       temp_group_current++;
       temp_group_compareTo++;
   } // All 2D planes have been processed
   
     //The back layer
        for (int r=0;r<rows;r++) {
            for (int c=0;c<columns-1;c++) { 
              /* We create a new border, with the appropriate pìxels */
              new_frontier->group1=temp_group_current;
              new_frontier->group2=temp_group_compareTo;
            
              /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
               if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
                distance=temp_group_current->reliability+temp_group_compareTo->reliability;
              else
                 distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
              if (mask!=NULL) {
                 if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                      distance=FLT_MAX;
              }

              //distance*=distance;
              new_frontier->distance=distance;
              new_frontier->distance_backup=0;
              new_frontier->is_banned=0;
              
              new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
            
              new_frontier++;
			
              temp_group_current++;
              temp_group_compareTo++;
              
            }
            temp_group_current++;
            temp_group_compareTo++;

        }
   
   
   
   if (debug)
     printf("  Z frontiers...\n");
   
      /* The borders of each pixels with the one in the inmediate 2D plane */
      temp_derivatives=derivatives;
      temp_group_current=array; // the current pixel
      temp_group_compareTo=array+twodsize; // the pixel to the right
      for (int d=0;d<depth-1;d++) {
        for (int i=0;i<rows-1;i++) {
          for (int j=0;j<columns-1;j++) {
            /* Update the values of the groups at both side of the border. Since there are only two pixels separated
            by a border, each group will be composed of a single pixel */
            new_frontier->group1=temp_group_current;
            new_frontier->group2=temp_group_compareTo;
          
            /* Calculate the distance */
            distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            if (qualities==NULL)
                /* The distance is given by the addition of the squares of the second derivatives in all directions */
              distance=temp_group_current->reliability+temp_group_compareTo->reliability;
            else
              distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
            if (mask!=NULL) {
              if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                  distance=FLT_MAX;
            }
            //distance*=distance;
            new_frontier->distance=distance;
            new_frontier->distance_backup=0;
            new_frontier->is_banned=0;
            
            new_frontier->steps=temp_derivatives->in_depth_steps;
          
            /* Insert the border in the heap */
            new_frontier++;
            temp_group_current++;
            temp_group_compareTo++;
            temp_derivatives++;
         } // a row has been processed except for the last column
         /* We create a new border, with the appropriate pixels */
         new_frontier->group1=temp_group_current;
         new_frontier->group2=temp_group_compareTo;
         
         /* The distance is given by the addition of the squares of the second derivatives in all directions */
         // The vertical steps of the last column were not calculated. They need to be calculated now
         distance=temp_group_current->reliability+temp_group_compareTo->reliability;
         if (qualities==NULL)
            /* The distance is given by the addition of the squares of the second derivatives in all directions */
           distance=temp_group_current->reliability+temp_group_compareTo->reliability;
         else
           distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
         if (mask!=NULL) {
           if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
               distance=FLT_MAX;
         }
        //distance*=distance;
         new_frontier->distance=distance;
         new_frontier->distance_backup=0;
         new_frontier->is_banned=0;
         
         new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
         
         /* insert the border in the heap */
         new_frontier++;
        
         temp_group_current++;
         temp_group_compareTo++;
         
         //temp_group_current++;
         //temp_group_compareTo++;
       } // a 2D plane has been processed except for the last row
       // There are no pre-calculated distances for the last row
       for (j=0;j<columns;j++) {
          /* We update the pointers with the groups, which are both composed of a single pixel */
          new_frontier->group1=temp_group_current;
          new_frontier->group2=temp_group_compareTo;
        
          /* Update the distance between the pixels at both side of the border */
          distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          if (qualities==NULL)
             /* The distance is given by the addition of the squares of the second derivatives in all directions */
             distance=temp_group_current->reliability+temp_group_compareTo->reliability;
          else
             distance=exp(-pow(*(qualities+(temp_group_current-array)),2)-pow(*(qualities+(temp_group_compareTo-array)),2));
              
          if (mask!=NULL) {
             if ((mask+(temp_group_current-array)==0) || (mask+(temp_group_compareTo-array)==0))
                 distance=FLT_MAX;
          }
         //distance*=distance;
          new_frontier->distance=distance;
          new_frontier->distance_backup=0;
          new_frontier->is_banned=0;
          
          new_frontier->steps=calculate_steps(temp_group_current->value,temp_group_compareTo->value);
        
          /* Insert the new border in the heap*/
          new_frontier++;
          temp_group_current++;
          temp_group_compareTo++;
          //temp_derivatives++;
       } // Now the entire 2D plane has been processed
       // Update pointers
       //temp_group_current++;
       //temp_group_compareTo++;
   } // All 2D planes have been processed
   
   if (debug) {
     printf("A total of %d frontiers have been created\n",new_frontier-frontier_array);
   }

}



/*************************************************************************************/
/*             FUNCTIONS RELATED TO SORT OPERATION ON THE BORDERS                    */  
/*************************************************************************************/
// It mixes two already sorted sections of the array
// Indexes of the first section starts at index1 and contains size1 elements
// Indexes of the second section starts at index1 and contains size1 elements
// frontierList is the array which contains the frontiers which are to be sorted 
// result is a temporary array that will act as an intermediate variable
void mixTogether(int *index1, int *index2, int size1,int size2, int* result)
{
  int counter1=0;
  int counter2=0;
  int *pointer_temp=index1;

  int *follower=result; // Where the final ordering will be stored
  //fprintf(f,"%d %d %d %d\n",*index1,*index2,size1,size2);
  //fflush(f);

  while ((counter1<size1) && (counter2<size2))
  // While we do not reach the end of one section
  // we chose the highest and add it to the temporal array follower, used to store the result
  {
	    if ((frontier_array[*(index1+counter1)].distance <= frontier_array[*(index2+counter2)].distance)) {
		// Depending on which element is the lowest
			// if the lowest is at the first portion
            *follower=*(index1+counter1);
			counter1++;
		} else {
          // If it is at the second portion     
          *follower=*(index2+counter2);
		  counter2++;
		}
		follower++; // increase the position of the array where we store the result
  }

  if (counter1==size1) { // Once we have finished the elements in one section we can just
                       // copy the rest as a block 
	  // We have terminated the first portion. We copy the rest of the second
	  // We copy it directly to the result
	  // First we copy the sorted section in follower
      memcpy(index1,result,sizeof(int)*(size1+counter2)); // copy as a block
	  // Then we copy the remaining of the second portion
      memcpy(index1+size1+counter2,(index2+counter2),sizeof(int)*(size2-counter2));
  } else {
	  // We have terminated the second portion. We copy the rest of the first
	  // In this case we have to it differently. Otherwise we would overwrite
	  // a part of the final result
	  memcpy(index1+size2+counter1,(index1+counter1),sizeof(int)*(size1-counter1));
	  // First we copy the sorted section in follower
      memcpy(index1,result,sizeof(int)*(size2+counter1)); // copy 
	  // Then we copy the remaining of the first portion
    }

  

}

void sort(int *index, int size, int* temp){
     if (size>1) {
		   int packetSize=size/2;
		   sort(index,packetSize,temp); // sort one half
		   sort((index+packetSize),size-packetSize,temp); // sort the other half
		   mixTogether(index,index+packetSize,packetSize,size-packetSize,temp); // merge
	}
}


void quickSort() {
	int *index=(int*)calloc(number_of_frontiers,sizeof(int)); // Array of indexes to work with positions
	int *temp=(int*)calloc(number_of_frontiers,sizeof(int)); // An array to support the merging operation
	         // It avoids a local allocation operation each time
    frontier* result_temp=(frontier*)calloc(number_of_frontiers, sizeof(frontier)); // to temporarily store the result
	frontier* current_result;
	int *current_index=index;
    int i;
	
    for (i=0;i<number_of_frontiers;i++) { // Initialize an array of indexes to work with the positions
      *current_index=i;
      current_index++;
    }

	sort(index,number_of_frontiers,temp);
	free(temp);       // temp is not further required

	
	current_result=result_temp;
	current_index=index;

    for (i=0;i<number_of_frontiers;i++){
	  *current_result=frontier_array[(*current_index)];
	  current_result++;
	  current_index++;
    }
    // Copy the data back to the original array
	memcpy(frontier_array,result_temp,number_of_frontiers*sizeof(frontier));

	free(index);
	free(result_temp);
}

int isNumber(char* string){
    int isNumber=1;
    while (isNumber==1 && *string!='\0') {
          if (*string<48 || *string>57)
            isNumber=0;
          string++;  
   }
   return isNumber;
}
    
              
/*************************************************************************************/
/*************************************************************************************/
/*                                                                                   */  
/*                               MAIN PROGRAM                                        */  
/*                                                                                   */  
/*************************************************************************************/
/*************************************************************************************/

int main(int argc, char *argv[])
{
   
   /*  variable declaration */

   group* group1;
   group* group2;
   int increment;
   frontier* current_frontier;
    
   clock_t start;
   clock_t finish;
   float duration;
   char* fileName;
   char* outputFileName;
   char* maskName;
   char* qualityName;
   int binaryFile;
   int number_of_borders_banned=0;
   int parameter_number=1;
   int parameter_index=1;
   int hasMask=0;
   int hasQualityMap=0;
   int error;
   char *parameter;
   mask=NULL;
   qualities=NULL;
   
   
   binaryFile=1;
   error=0;
   while (!error && argv[parameter_index]!=NULL) //parameter_number<5)
   {
         parameter=argv[parameter_index];
         //printf("parameter:  %s \n",parameter);       
         if (strcmp("-t",parameter)==0) {
           //printf("-t\n");
           binaryFile=0;
           parameter_index++;
         } else
         if (strcmp("-m",parameter)==0) {
           //printf("-m\n");
           parameter_index++;                   
           parameter=argv[parameter_index];
           if (parameter==NULL)
             error=true;
           else {
             maskName=parameter;
             hasMask=1;
           }
           parameter_index++;
         } else
         if (strcmp("-q",parameter)==0){
           //printf("-q\n");
           parameter_index++;                   
           parameter=argv[parameter_index];
           if (parameter==NULL)
             error=true;
           else {
             qualityName=parameter;
             hasQualityMap=1;
           }
           parameter_index++;
         } else 
          {
             //printf("par\n");
             switch (parameter_number) {
                    case 1: fileName=argv[parameter_index]; break;
                    case 2: outputFileName=argv[parameter_index]; break;
                    case 3: if (!isNumber(argv[parameter_index])) 
                               error = true;
                            else
                              rows=atoi(argv[parameter_index]); break;
                    case 4: if (!isNumber(argv[parameter_index])) 
                               error = true;
                            else
                              columns=atoi(argv[parameter_index]); break;
                    case 5: if (!isNumber(argv[parameter_index])) 
                               error = true;
                            else
                              depth=atoi(argv[parameter_index]); break;
                    default: error=true;
             }
             parameter_index++;
             parameter_number++;
          }
             
   }

   // Check that the parameters are correct
   if (error) 
   {
	   printf("______________________________________________________________________________\n");
       printf("\nThe program syntax is:\n                  unwrap [-t] [-m mask_file] [-q quality_map_file} input_file_name ouput_file_name n_rows n_columns depth\n");
       printf("______________________________________________________________________________\n");
       printf("\ninput_file_name - file which contains the data to be unwrapped\n");
       printf("output_file_name - output file to write the result to\n");
       printf("n_rows - number of rows in the image\n");
       printf("n_columns - number of columns in the image\n");
       printf("n_depth - depth of the image\n");
       printf("\n\nOPTIONS\n");
       printf("-------\n");
       printf("use unwrap -t to use plain text files\n");
       printf("use -m mask_file to specify a mask file (unsigned chars with 0 values at the mask)\n");
       printf("use -q quality_map_file to specify a quality map (otherwise second differences are used)\n");
       printf("______________________________________________________________________________\n");
       printf("\nALGORITHM DEVELOPED IN A COLLABORATION BETWEEN:\n");
       printf("GERI (GENERAL ENGINEERING RESEARCH INSTITUTE) - LIVERPOOL JOHN MOORES UNIVERSITY - UK\n");
       printf("UNIVERSITY OF VALENCIA (COMPUTING DEPARTMENT) - SPAIN\n");
       printf("______________________________________________________________________________\n");
       printf("\nIMPORTANT NOTE: If you use this algorithm, we would like to be informed. Please\n");
       printf("write an e-mail to miguel.arevalillo@uv.es or m.a.gdeisat@ljmu.ac.uk\n");
       printf("We take absolute no responsability of any damage derived from the use of this code\n");
       printf("______________________________________________________________________________\n");
     //printf("Number of parameters: %d\n",parameter_index); 
     //printf("Number of parameters: %d\n",parameter_number); 

	   exit(0);
   }
    
   number_of_pixels=rows*columns*depth;
   if (debug) {
              
     printf("Unwrapping %d x %d x %d 3D image in file %s\n",rows,columns,depth,fileName );
     printf("Total number of pixels: %d\n",number_of_pixels); 
     if (hasQualityMap)
        printf("Quality map file: %s\n",qualityName); 
     if (hasMask)
        printf("Mask file: %s\n",maskName); 
     if (binaryFile==1)
        printf("Files are binary\n");
     else
        printf("Files are in text format\n");
     
     //printf("Number of parameters: %d\n",parameter_index); 
     //printf("Number of parameters: %d\n",parameter_number); 
     
   }
   
   if (hasMask) {
     if (debug)
       printf("Loading mask file (%s)... ",maskName);
     mask=(unsigned char*)calloc(number_of_pixels,sizeof(float));
     if (binaryFile)
       load_mask_binary(fileName,mask);
     else
       load_mask_plain_text(fileName,mask);
   }
   
   if (hasQualityMap) {
     if (debug)
       printf("Loading quality map file (%s)... ",qualityName);
     qualities=(float*)calloc(number_of_pixels,sizeof(float));
     
     if (binaryFile)
       load_file_binary(fileName,qualities);
     else
       load_file_plain_text(fileName,qualities);
   }
                      

   values_array=(float*)calloc(number_of_pixels,sizeof(float));
   /*  Load the file  */
   if (debug)
     printf("Loading input file... ");
   if (binaryFile)
     load_file_binary(fileName,values_array);
   else
     load_file_plain_text(fileName,values_array);
   
   printf(" done\n");
   /* Start the clock */
   if (debug)
     printf("Allocating memory... ");
   start = clock(); 
   
   /*  Reserve memory for the pixels  */
   array=(group*)calloc(number_of_pixels,sizeof(group));
   convert_values_to_struct(); // Transfer the values to the struct. This function has been implemented
                               // separate to allow for execution time measurement
   free(values_array);
   
   twodsize=rows*columns;
   r_twodsize=(rows-1)*(columns-1);
   number_of_frontiers= (columns-1)*rows*depth+(rows-1)*columns*depth+rows*columns*(depth-1);
   
   
   
   /*  Reserve memory for the 2 x 2 loops  */
   residues=(residue_values*)calloc(number_of_pixels,sizeof(residue_values));
   
   frontier_array=(frontier*)calloc(number_of_frontiers, sizeof(frontier));

   hfrontiers=frontier_array+(rows-1)*columns*depth; //   x | x 
   dfrontiers=hfrontiers+rows*(columns-1)*depth;
   
   
   if (debug)
     printf("done\n");
   
   if (debug)
     printf("Image size: %d\n",number_of_pixels);
   
   if (debug)
     printf("Preparing data structures...");
   
      
   /*  Reserve the vector to store the distances in all directions */
   derivatives=(distances*)calloc((columns-1)*(rows-1)*(depth-1),sizeof(distances));
   if (debug)
     printf(" done\n");
   
   /****** Prepare the data structure:  *****/
   
   /*  Calculate pixel reliability    */
   calculate_pixel_reliability();

  
   /*  Build the borders */
   if (debug)
     printf("Registering all possible borders...\n");

   if ((qualities!=NULL) || (mask!=NULL)) {
     build_frontiers();
   } else 
      build_frontiers(mask,qualities);

   /*  Free the space used to store the distances */
   free(derivatives);
   
   if (debug)
     printf("Processing residue loops and discarding borders...");
   build_and_process_all_loops();
   if (debug)
     printf(" done\n");

  
   printf("sorting by quality...");
   quickSort();
   if (debug)
     printf("done\n");

      
   
   if (debug)
     printf("Starting data processing... \n");
   /*  While there are elements in the heap */
   
   current_frontier=frontier_array;
   
   if (debug)
     printf("Processing a total of %d borders...",number_of_frontiers);
   
   for (long i=0;i<number_of_frontiers;i++)
   {
	  if (!current_frontier->is_banned){
      // If the frontier was banned we do not process it
        group1=getAndSortGroup(current_frontier->group1);
        group2=getAndSortGroup(current_frontier->group2);
	    if (group1!=group2) {
		  increment=current_frontier->steps;
          // increment stores the number of 2pi jumps which should be added to all elements
          // in group2 to achieve a local unwrapping. group2 is concatenated to group1 
          group2->current_increment=current_frontier->group1->current_increment-
			   current_frontier->group2->current_increment+increment;
		  group2->previous_group=group1;
          /* Join the pixels in both groups */
          group1->last_pixel->next_pixel=group2;
          group1->last_pixel=group2->last_pixel;
        } 
      } else
         number_of_borders_banned++;
      /* Free the memory of the border */
      current_frontier++;
    } // end of the while that controls the heap
   if (debug)
     printf(" done\n");

   if (debug)
     printf("A total of %d borders were banned by the algorithm\n",number_of_borders_banned);
    
    free(frontier_array); 
	
   /* Calculate the final values */
    if (debug)
       printf("Calculating unwrapped final values...");
   	calculate_final_pixel_values();
    if (debug)
      printf(" done\n");
     
	/* Stop the clock */
	finish = clock();
	duration = (float)(finish - start) / CLOCKS_PER_SEC;
	
	/* Print the processing time */
    if (debug)
      printf( "Data processed in %2.1f seconds.\nNote that this time includes the display of messages\n", duration );
	
    if (debug)
       printf("Saving results in file...");
    /* Store the results */
    if (binaryFile)
       save_file_binary(outputFileName,array); 
    else
       save_file_plain_text(outputFileName,array); 
    //save_file_plain_text(outputFileName,array); 
    if (debug)
      printf(" done\n");
    if (debug)
       printf("The unwrapped data has been stored in %s\n",outputFileName);
    
    /* Free the space */
    free(array);
};
