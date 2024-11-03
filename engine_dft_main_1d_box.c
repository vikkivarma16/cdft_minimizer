#include <stdio.h>
#include <stdlib.h>
#include "cjson/cJSON.h"
#include <string.h>  // Fixed the include directive for <string.h>
#include <fftw3.h>
#include <math.h>


///##########################################################################///

/// variable declaration .......................................................

int i, j, k, l, m, n;

double bdm, dummy, temp;

const double piee = acos(-1);



//###########################################################################///
/// Data reading sections....
////////////////////////////////////////////////////////////////////////////////
// Defining the box properties

typedef struct {
    int dimension;
    char confinement[10];
} space_data;

typedef struct {
    double x;
    double y;
    double z;
} box_length;

typedef struct {
    int nx;
    int ny;
    int nz;
} box_points;

typedef struct {
    box_length length;
    box_points points;
} box_data;
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Define the thermodynamic properties
typedef struct {
    double temperature;
    double rho;
    int iteration_max;
} thermodynamic_data;
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Define structure for interaction properties
typedef struct {
    char id[10];
    double rho_frac;
} species_data;
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Define structure for categories of interactions
typedef struct {
    char id[10];
    char type[10];  // Assuming type is a short string (e.g., "lj", "gs")
    double sigma;
    double cutoff;
    double epsilon;
} interaction_pair;

// Define structure for all interactions
typedef struct {
    interaction_pair primary;
    interaction_pair secondary;
    interaction_pair tertiary;
} interaction_data;
////////////////////////////////////////////////////////////////////////////////

///##########################################################################///










///##########################################################################///

////////////////////////////////////////////////////////////////////////////////
// Function to parse species data
void parse_data_species(cJSON *species_properties, species_data *species) {
    cJSON *element_species;
    int i = 0;  // Initialize i to avoid potential issues
    cJSON_ArrayForEach(element_species, species_properties) {
        strcpy(species[i].id, element_species->string);
        species[i].rho_frac = cJSON_GetObjectItem(element_species, "rho_frac")->valuedouble;
        i = i + 1;
    }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Function to parse interaction data
void parse_interactions_data(cJSON *interactions_properties, interaction_data *interactions, int n_interactions) {
    cJSON *primary_interactions = cJSON_GetObjectItem(interactions_properties, "primary");

    int i, flag;
    cJSON *pair;
    
    for (i = 0; i < n_interactions; i++) {
        flag = 0;
        cJSON_ArrayForEach(pair, primary_interactions) 
        {
            
            if (strcmp(pair->string, interactions[i].primary.id) == 0) {
                flag = 1;
                ///printf("Primary is running dear %s \n\n", pair->string);
                cJSON *temp_data = cJSON_GetObjectItem(primary_interactions, pair->string);
                strcpy(interactions[i].primary.type, cJSON_GetObjectItem(temp_data, "type")->valuestring);
                interactions[i].primary.sigma = cJSON_GetObjectItem(temp_data, "sigma")->valuedouble;
                interactions[i].primary.cutoff = cJSON_GetObjectItem(temp_data, "cutoff")->valuedouble;
                interactions[i].primary.epsilon = cJSON_GetObjectItem(temp_data, "epsilon")->valuedouble;
            }
            
        }
        
        if (flag == 0) 
        {
              strcpy(interactions[i].primary.type, "NA");
              interactions[i].primary.sigma = 0.0;
              interactions[i].primary.cutoff = 0.0;
              interactions[i].primary.epsilon = 0.0;
        }
    }
    
    cJSON *secondary_interactions = cJSON_GetObjectItem(interactions_properties, "secondary");
    
    
    for (i = 0; i < n_interactions; i++) {
        flag = 0;
        cJSON_ArrayForEach(pair, secondary_interactions) 
        {
            
            if (strcmp(pair->string, interactions[i].primary.id) == 0) {
                flag = 1;
                ///printf("Secondary is running dear %s \n\n", pair->string);
                cJSON *temp_data = cJSON_GetObjectItem(secondary_interactions, pair->string);
                strcpy(interactions[i].secondary.type, cJSON_GetObjectItem(temp_data, "type")->valuestring);
                interactions[i].secondary.sigma = cJSON_GetObjectItem(temp_data, "sigma")->valuedouble;
                interactions[i].secondary.cutoff = cJSON_GetObjectItem(temp_data, "cutoff")->valuedouble;
                interactions[i].secondary.epsilon = cJSON_GetObjectItem(temp_data, "epsilon")->valuedouble;
            }
            
        }
        
        if (flag == 0) 
        {
              strcpy(interactions[i].secondary.type, "NA");
              interactions[i].secondary.sigma = 0.0;
              interactions[i].secondary.cutoff = 0.0;
              interactions[i].secondary.epsilon = 0.0;
        }
    }
    
    cJSON *tertiary_interactions = cJSON_GetObjectItem(interactions_properties, "tertiary");
    for (i = 0; i < n_interactions; i++) {
        flag = 0;
        cJSON_ArrayForEach(pair, tertiary_interactions) 
        {
            
            if (strcmp(pair->string, interactions[i].primary.id) == 0) {
                flag = 1;
                ///printf("tertiary is running dear %s \n\n", pair->string);
                cJSON *temp_data = cJSON_GetObjectItem(tertiary_interactions, pair->string);
                strcpy(interactions[i].tertiary.type, cJSON_GetObjectItem(temp_data, "type")->valuestring);
                interactions[i].tertiary.sigma = cJSON_GetObjectItem(temp_data, "sigma")->valuedouble;
                interactions[i].tertiary.cutoff = cJSON_GetObjectItem(temp_data, "cutoff")->valuedouble;
                interactions[i].tertiary.epsilon = cJSON_GetObjectItem(temp_data, "epsilon")->valuedouble;
            }
            
        }
        
        if (flag == 0) 
        {
              strcpy(interactions[i].tertiary.type, "NA");
              interactions[i].tertiary.sigma = 0.0;
              interactions[i].tertiary.cutoff = 0.0;
              interactions[i].tertiary.epsilon = 0.0;
        }
    }
    
    
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Function to parse box length
void parse_data_box_length(cJSON *box_prop_length, box_length *b_length) {
    cJSON *bx = cJSON_GetArrayItem(box_prop_length, 0);
    cJSON *by = cJSON_GetArrayItem(box_prop_length, 1);
    cJSON *bz = cJSON_GetArrayItem(box_prop_length, 2);

    if (bx != NULL) b_length->x = bx->valuedouble;
    if (by != NULL) b_length->y = by->valuedouble;
    if (bz != NULL) b_length->z = bz->valuedouble;
}

// Function to parse box points
void parse_data_box_points(cJSON *box_prop_points, box_points *b_points) {
    cJSON *bx = cJSON_GetArrayItem(box_prop_points, 0);
    cJSON *by = cJSON_GetArrayItem(box_prop_points, 1);
    cJSON *bz = cJSON_GetArrayItem(box_prop_points, 2);

    if (bx != NULL) b_points->nx = bx->valueint;
    if (by != NULL) b_points->ny = by->valueint;
    if (bz != NULL) b_points->nz = bz->valueint;
}
////////////////////////////////////////////////////////////////////////////////

///##########################################################################///








///##########################################################################///

////////////////////////////////////////////////////////////////////////////////
///... Calculator sections ..................................................///


void calculate_rho_fmt_1d(fftw_complex ***k_data_fmt_1d, fftw_complex ***r_data_fmt_1d, fftw_complex ***kernel_data, fftw_plan *fftw_plan_forward_rho, fftw_plan **fftw_plan_backward_rho, fftw_plan **fftw_plan_forward_dphi, fftw_plan **fftw_plan_backward_dphi, int N_kx, int n_species, int *fmt_flag){
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
           fftw_execute(fftw_plan_forward_rho[i]);
      }
      //////////////////////////////////////////////////////////////////////////
     
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
            k=13;      
            for (j = 0 ; j < 4; j++)
            {
                
                for (l = 0; l < N_kx; l++)
                {
                    k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][0][l][0] * kernel_data[i][j][l][0];
                    k_data_fmt_1d[i][k][l][1] = k_data_fmt_1d[i][0][l][1] * kernel_data[i][j][l][0];
                }
                k=k+1;
            }
            
            for (l = 0; l < N_kx; l++)
            {
                k_data_fmt_1d[i][k][l][1] = -k_data_fmt_1d[i][0][l][0] * kernel_data[i][4][l][0];
                k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][0][l][1] * kernel_data[i][4][l][0];
            }
            
            
            k=k+1;
            for (l = 0; l < N_kx; l++)
            {
                k_data_fmt_1d[i][k][l][1] = -k_data_fmt_1d[i][0][l][0] * kernel_data[i][5][l][0];
                k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][0][l][1] * kernel_data[i][5][l][0];
            }
           
            
      }
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
           for (j=0; j<6; j++)
           {
                fftw_execute(fftw_plan_backward_rho[i][j]);
           }
      }
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      double total_n_alpha[6][2];
      for (l = 0; l < N_kx; l++)
      { 
          
          for (i = 0; i<6; i++)
          {
              total_n_alpha[i][0]=0.0;
              total_n_alpha[i][1]=0.0;
               
              for(k=0; k<n_species; k++)
              {
                  total_n_alpha[i][0]=total_n_alpha[i][0]+r_data_fmt_1d[k][i+1][l][0];
                  total_n_alpha[i][1]=total_n_alpha[i][1]+r_data_fmt_1d[k][i+1][l][1];
              }        
          }
          for (k = 0; k < n_species; k++)
          {
              r_data_fmt_1d[k][7][l][0] = -log(1.0-total_n_alpha[3][0]);
              r_data_fmt_1d[k][7][l][1] = 0.0;      
                        
              r_data_fmt_1d[k][8][l][0] = total_n_alpha[2][0]/(1.0-total_n_alpha[3][0]);
              r_data_fmt_1d[k][8][l][1] = 0.0;             
              
              r_data_fmt_1d[k][9][l][0] = total_n_alpha[1][0]/(1.0-total_n_alpha[3][0]) + (1.0/(8.0*piee))*(total_n_alpha[2][0]*total_n_alpha[2][0]- total_n_alpha[4][0]*total_n_alpha[4][0]) /((1.0-total_n_alpha[3][0]) * (1.0-total_n_alpha[3][0]));
              r_data_fmt_1d[k][9][l][1] = 0.0;    
              
              r_data_fmt_1d[k][10][l][0] = total_n_alpha[0][0]/(1.0-total_n_alpha[3][0])  +  (total_n_alpha[1][0] * total_n_alpha[2][0]-total_n_alpha[4][0]*total_n_alpha[5][0])/((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0])) ;
              r_data_fmt_1d[k][10][l][0] = r_data_fmt_1d[k][10][l][0] + (1.0/(12.0*piee))*(total_n_alpha[2][0]*total_n_alpha[2][0]*total_n_alpha[2][0]-3.0*total_n_alpha[5][0]*total_n_alpha[5][0])/ ((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0]));
              r_data_fmt_1d[k][10][l][1] = 0.0;            
              
              r_data_fmt_1d[k][11][l][0] = total_n_alpha[5][0]/(1.0-total_n_alpha[3][0]);
              r_data_fmt_1d[k][11][l][1] = 0.0;    
              
              r_data_fmt_1d[k][12][l][0] = total_n_alpha[4][0]/(1.0-total_n_alpha[3][0]) + 0.25*(total_n_alpha[2][0]*total_n_alpha[5][0])/((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0]));
              r_data_fmt_1d[k][12][l][1] = 0.0;    
                                          
          }
            
      }
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
          for (j = 0; j < 6; j++)
          {
              fftw_execute(fftw_plan_forward_dphi[i][j]);
          }
      }     
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
            k=13;      
            for (j = 0 ; j < 4; j++)
            {
                
                for (l = 0; l < N_kx; l++)
                {
                    k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][j+7][l][0] * kernel_data[i][j][l][0];
                    k_data_fmt_1d[i][k][l][1] = k_data_fmt_1d[i][j+7][l][1] * kernel_data[i][j][l][0];
                }
                k=k+1;
            }
            

            for (l = 0; l < N_kx; l++)
            {
                k_data_fmt_1d[i][k][l][1] = -k_data_fmt_1d[i][j+7][l][0] * kernel_data[i][4][l][0];
                k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][j+7][l][1] * kernel_data[i][4][l][0];
            }
            
            
            k=k+1;
            j=j+1;
            for (l = 0; l < N_kx; l++)
            {
                k_data_fmt_1d[i][k][l][1] = -k_data_fmt_1d[i][j+7][l][0] * kernel_data[i][5][l][0];
                k_data_fmt_1d[i][k][l][0] = k_data_fmt_1d[i][j+7][l][1] * kernel_data[i][5][l][0];
            } 
      }
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////
      for (i = 0; i < n_species; i++)
      {
          for (j = 0; j < 6; j++)
          {
              fftw_execute(fftw_plan_backward_dphi[i][j]);
          }
      }    
      //////////////////////////////////////////////////////////////////////////
      
      //////////////////////////////////////////////////////////////////////////

}








////////////////////////////////////////////////////////////////////////////////





















///MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM///


int main() {



    ///######################################################################///
    
    ////////////////////////////////////////////////////////////////////////////
    FILE *file_space = fopen("input_space_properties.json", "r");
    if (file_space == NULL) {
        perror("Could not open the file");
        return EXIT_FAILURE;
    }
    fseek(file_space, 0, SEEK_END);
    long length_space = ftell(file_space);
    fseek(file_space, 0, SEEK_SET);

    char *data_space = malloc(length_space + 1);
    fread(data_space, 1, length_space, file_space);
    data_space[length_space] = '\0';

    cJSON *json_space = cJSON_Parse(data_space);
    cJSON *space_properties = cJSON_GetObjectItem(json_space, "space_properties");

    // to store the space data in the structure format
    space_data space;
    space.dimension = cJSON_GetObjectItem(space_properties, "dimension")->valueint;
    strcpy(space.confinement, cJSON_GetObjectItem(space_properties, "confinement")->valuestring);

    fclose(file_space);
    cJSON_Delete(json_space);
    free(data_space);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    FILE *file_box = fopen("input_box_properties.json", "r");
    if (file_box == NULL) {
        perror("There is no file named input_box_properties.json");
        return EXIT_FAILURE;
    }
    fseek(file_box, 0, SEEK_END);
    long length_box = ftell(file_box);
    fseek(file_box, 0, SEEK_SET);

    char *data_box = malloc(length_box + 1);
    fread(data_box, 1, length_box, file_box);
    data_box[length_box] = '\0';
    fclose(file_box);

    cJSON *json_box = cJSON_Parse(data_box);
    if (json_box == NULL) {
        perror("JSON parsing failed. Please check for inconsistencies.");
        return EXIT_FAILURE;
    }

    box_data box = { .length = {0.0, 0.0, 0.0}, .points = {0, 0, 0} };

    cJSON *box_prop = cJSON_GetObjectItem(json_box, "box_properties");
    parse_data_box_length(cJSON_GetObjectItem(box_prop, "box_length"), &box.length);
    parse_data_box_points(cJSON_GetObjectItem(box_prop, "box_points"), &box.points);

    cJSON_Delete(json_box);
    free(data_box);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    FILE *file_thermodynamic_data = fopen("input_thermodynamic_properties.json", "r");
    if (file_thermodynamic_data == NULL) {
        perror("Could not open the thermodynamic file");
        return EXIT_FAILURE;
    }

    fseek(file_thermodynamic_data, 0, SEEK_END);
    long length_ftd = ftell(file_thermodynamic_data);
    fseek(file_thermodynamic_data, 0, SEEK_SET);

    char *data_thermodynamic = malloc(length_ftd + 1);
    fread(data_thermodynamic, 1, length_ftd, file_thermodynamic_data);
    data_thermodynamic[length_ftd] = '\0';
    fclose(file_thermodynamic_data);

    cJSON *json_thermodynamic = cJSON_Parse(data_thermodynamic);
    cJSON *tp = cJSON_GetObjectItem(json_thermodynamic, "thermodynamic_properties");

    thermodynamic_data thermo;
    thermo.temperature = cJSON_GetObjectItem(tp, "temperature")->valuedouble;
    thermo.rho = cJSON_GetObjectItem(tp, "rho")->valuedouble;
    thermo.iteration_max = cJSON_GetObjectItem(tp, "iteration_max")->valueint;

    cJSON_Delete(json_thermodynamic);
    free(data_thermodynamic);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    FILE *file_species = fopen("input_species_properties.json", "r");
    if (file_species == NULL) {
        perror("Could not open species file");
        return EXIT_FAILURE;
    }

    fseek(file_species, 0, SEEK_END);
    long length_species = ftell(file_species);
    fseek(file_species, 0, SEEK_SET);

    char *data_species = malloc(length_species + 1);
    fread(data_species, 1, length_species, file_species);
    data_species[length_species] = '\0';
    fclose(file_species);

    cJSON *json_species = cJSON_Parse(data_species);
    cJSON *species_properties = cJSON_GetObjectItem(json_species, "species");

    int n_species = cJSON_GetArraySize(species_properties);
    species_data *species = malloc(n_species * sizeof(species_data));
    parse_data_species(species_properties, species);

    cJSON_Delete(json_species);
    free(data_species);
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    FILE *file_interactions = fopen("input_interactions_properties.json", "r");
    if (file_interactions == NULL) {
        perror("Could not open interactions file");
        return EXIT_FAILURE;
    }

    fseek(file_interactions, 0, SEEK_END);
    long length_interactions = ftell(file_interactions);
    fseek(file_interactions, 0, SEEK_SET);

    char *data_interactions = malloc(length_interactions + 1);
    fread(data_interactions, 1, length_interactions, file_interactions);
    data_interactions[length_interactions] = '\0';
    fclose(file_interactions);

    cJSON *json_interactions = cJSON_Parse(data_interactions);
    cJSON *interactions_properties = cJSON_GetObjectItem(json_interactions, "interactions");

    int n_interactions; 
    n_interactions = (int) (n_species*(n_species-1)*0.5) + n_species;
    interaction_data *interactions = malloc(n_interactions * sizeof(interaction_data));
    
    
    char part1[10], part2[10], total_char[10];
    
    int k = 0;
    for (int i = 0; i < n_species; i++) {
        for (int j = i; j < n_species; j++) {
            strcpy(part1, "");
            strcpy(part2, "");
            strcpy(total_char, "");

            // Copy species IDs into part1 and part2
            strcpy(part1, species[i].id);
            strcpy(part2, species[j].id);

            // Concatenate part1 and part2 into total
            strcpy(total_char, part1);
            strcat(total_char, part2);

            // Copy total into interactions[k].primary.id
            strcpy(interactions[k].primary.id, total_char);
            k++;
        }
    }
    
    
    k = 0;
    for (int i = 0; i < n_species; i++) {
        for (int j = i; j < n_species; j++) {
            strcpy(part1, "");
            strcpy(part2, "");
            strcpy(total_char, "");

            // Copy species IDs into part1 and part2
            strcpy(part1, species[i].id);
            strcpy(part2, species[j].id);

            // Concatenate part1 and part2 into total
            strcpy(total_char, part1);
            strcat(total_char, part2);

            // Copy total into interactions[k].primary.id
            strcpy(interactions[k].secondary.id, total_char);
            k++;
        }
    }

     k = 0;
    for (int i = 0; i < n_species; i++) {
        for (int j = i; j < n_species; j++) {
            strcpy(part1, "");
            strcpy(part2, "");
            strcpy(total_char, "");

            // Copy species IDs into part1 and part2
            strcpy(part1, species[i].id);
            strcpy(part2, species[j].id);

            // Concatenate part1 and part2 into total
            strcpy(total_char, part1);
            strcat(total_char, part2);

            // Copy total into interactions[k].primary.id
            strcpy(interactions[k].tertiary.id, total_char);
            k++;
        }
    }

    
    parse_interactions_data(interactions_properties, interactions, n_interactions);

    cJSON_Delete(json_interactions);
    free(data_interactions);
    ////////////////////////////////////////////////////////////////////////////
    
    ///######################################################################///
    
    
    
    
    
    
    
    
    
    
    
    ///######################################################################///
    
    ////////////////////////////////////////////////////////////////////////////
    
    
    // Print space data
    printf("\n\n ... congratulations for parsing all the data ... your system have been fed with the following information ... \n\n\n");
    printf("\n\n\n... space data: ...\n");
    printf("Dimension: %d\n", space.dimension);
    printf("Confinement: %s\n\n\n", space.confinement);

  
    // Print box data
    printf("... box data: ...\n");
    printf("Length: %lf, %lf, %lf\n", box.length.x, box.length.y, box.length.z);
    printf("Points: %d, %d, %d\n\n", box.points.nx, box.points.ny, box.points.nz);
    printf("\n\n\n");


     // Print thermodynamic data
    printf("... thermodynamic data: ... \n");
    printf("Temperature: %lf\n", thermo.temperature);
    printf("Density (rho): %lf\n", thermo.rho);
    printf("Max Iterations: %d\n", thermo.iteration_max);
    printf("\n\n\n");

    

    // Print species data
    printf("... species data: ...\n");
    for (int i = 0; i < n_species; i++) {
        printf("ID: %s, Density Fraction: %lf\n", species[i].id, species[i].rho_frac);
    }
    printf("\n\n\n");

    // Print interaction data
    printf("... interaction Data for the primary interactions: ...\n");
    for (int i = 0; i < n_interactions; i++) {
        printf("Interaction Pair ID: %s\n", interactions[i].primary.id);
        printf("  Type: %s\n", interactions[i].primary.type);
        printf("  Sigma: %lf\n", interactions[i].primary.sigma);
        printf("  Cutoff: %lf\n", interactions[i].primary.cutoff);
        printf("  Epsilon: %lf\n", interactions[i].primary.epsilon);
        printf("\n");
    }
    
    printf("\n\n");
    
    // Print interaction data
    printf("... interaction Data for the secondary interactions: ...\n");
    for (int i = 0; i < n_interactions; i++) {
        printf("Interaction Pair ID: %s\n", interactions[i].secondary.id);
        printf("  Type: %s\n", interactions[i].secondary.type);
        printf("  Sigma: %lf\n", interactions[i].secondary.sigma);
        printf("  Cutoff: %lf\n", interactions[i].secondary.cutoff);
        printf("  Epsilon: %lf\n", interactions[i].secondary.epsilon);
        printf("\n");
    }
    
    
    printf("\n\n");
    
    // Print interaction data
    printf("... interaction Data for the tertiary interactions: ...\n");
    for (int i = 0; i < n_interactions; i++) {
        printf("Interaction Pair ID: %s\n", interactions[i].tertiary.id);
        printf("  Type: %s\n", interactions[i].tertiary.type);
        printf("  Sigma: %lf\n", interactions[i].tertiary.sigma);
        printf("  Cutoff: %lf\n", interactions[i].tertiary.cutoff);
        printf("  Epsilon: %lf\n", interactions[i].tertiary.epsilon);
        printf("\n");
    }
    

    printf("\n\n\n");
    ////////////////////////////////////////////////////////////////////////////
    
    ///######################################################################///
   
    
    
    
    
    
    
    ///######################################################################///
    
    ////////////////////////////////////////////////////////////////////////////
    /////Now the main code for the calculation of the dft minimization began////
    
    
    
    
    
    
    
    
    
    ///space.dimension =2  ;
    if (space.dimension == 1)
    {
    
        ///##################################################################///
        
        ////////////////////////////////////////////////////////////////////////
        /// variable assignment sections ... 
        /// total variables assigned on each points including six kinds of densities, six kinds of dphis and the six kinds of rhos + one kind of real space total rho....
        
        int N_kx;
        
        int N_x;
        
        N_kx = box.points.nx;
        N_x = N_kx;
        
        
        
        int t_variables;
        t_variables = 1+6+6+6+1;  ///total_rho + rho_fmt + dphi+ fmt_replica; used for the multiplication + r_space_value for 1d space ..... ..... ..... ..... ..... ..... ..... ..... ..... .....
        
        fftw_complex***    k_data_fmt_1d = malloc(n_species*sizeof(fftw_complex**));
        
        for (i = 0; i< n_species; i++)
        {
            k_data_fmt_1d[i] = malloc(t_variables * sizeof (fftw_complex*));
            
            for (j=0; j < t_variables; j++)
            {
                k_data_fmt_1d[i][j] = malloc(N_kx * sizeof (fftw_complex));
                
                for (k = 0; k < N_kx; k++) {
                    k_data_fmt_1d[i][j][k][0] = 0.0; // Real part
                    k_data_fmt_1d[i][j][k][1] = 0.0; // Imaginary part
                }
            }
        
        }
        
        
        fftw_complex***    r_data_fmt_1d = malloc(n_species*sizeof(fftw_complex**));
        
        for (i = 0; i< n_species; i++)
        {
            r_data_fmt_1d[i] = malloc(t_variables * sizeof (fftw_complex*));
            
            for (j=0; j < t_variables; j++)
            {
                r_data_fmt_1d[i][j] = malloc(N_x * sizeof (fftw_complex));
                
                for (int k = 0; k < N_x; k++) {
                   r_data_fmt_1d[i][j][k][0] = 0.0; // Real part
                   r_data_fmt_1d[i][j][k][1] = 0.0; // Imaginary part
                }
            }
        }
        
        
        
        int t_kernel_variables;
        t_kernel_variables = 6;
        
        fftw_complex***  kernel_data = malloc(n_species*sizeof(fftw_complex**));
        
         for (i = 0; i< n_species; i++)
        {
            kernel_data[i] = malloc(t_kernel_variables*sizeof(fftw_complex*));
            for (j=0; j < t_kernel_variables; j++)
            {
               kernel_data[i][j] = malloc(N_kx*sizeof(fftw_complex));
               for (int k = 0; k < N_kx; k++) {
                    kernel_data[i][j][k][0] = 0.0; // Real part
                    kernel_data[i][j][k][1] = 0.0; // Imaginary part
               }
            }
        }
        
        ////////////////////////////////////////////////////////////////////////
        
        ///##################################################################///
        
        
        
        
        
        
        
        
       ///###################################################################///
       
       /////////////////////////////////////////////////////////////////////////
       /// input data being fed for the fmt kernels ....
       
       char temp_char[10] = {0};
       strcpy(temp_char, "hc");
       int* flag_fmt = malloc(n_species*sizeof(int));
       
       for (i = 0; i<n_species; i++)
       {
          flag_fmt[i]=0;
       }
       
       for (i = 0; i<n_species; i++)
       {
            k=0;
            for (j= 0; j<i; j++)
            {
                k=k+(n_species-j);
            
            }
            if (strcmp(interactions[k].primary.type, temp_char) == 0)
            {
                char filename_fmt_kernel[255] = {0};
                sprintf(filename_fmt_kernel, "weight_FMT_k_space_%s.txt", species[i].id);
                FILE* file_kernel = fopen(filename_fmt_kernel, "r");
                
                for (j = 0; j <N_kx; j++)
                {
                      for (l = 0; l < 3; l++)
                      {
                          fscanf(file_kernel, "%lf", &bdm);
                      }
                      for (l = 0; l < t_kernel_variables; l++)
                      {
                          fscanf(file_kernel, "%lf", &kernel_data[i][l][j][0]);
                      }
                }
                
                fclose(file_kernel);
            }
       }
       
       
       
       for (i = 0; i<n_species; i++)
       {
            k=0;
            for (j= 0; j<i; j++)
            {
                k=k+(n_species-j);
            
            }
  
            if (strcmp(interactions[k].secondary.type, temp_char) == 0)
            {
                char filename_fmt_kernel[255] = {0};
                sprintf(filename_fmt_kernel, "weight_FMT_k_space_%s.txt", species[i].id);
                FILE* file_kernel = fopen(filename_fmt_kernel, "r");
                
                for (j = 0; j <N_kx; j++)
                {
                      for (l = 0; l < 3; l++)
                      {
                          fscanf(file_kernel, "%lf", &bdm);
                      }
                      for (l = 0; l < t_kernel_variables; l++)
                      {
                          fscanf(file_kernel, "%lf", &kernel_data[i][l][j][0]);
                      }
                }
                fclose(file_kernel);
            }
       }
       
      
       for (i = 0; i<n_species; i++)
       {
            k=0;
            for (j= 0; j<i; j++)
            {
                k=k+(n_species-j);
            
            }
            
            if (strcmp(interactions[k].tertiary.type, temp_char) == 0)
            {
                char filename_fmt_kernel[255] = {0};
                sprintf(filename_fmt_kernel, "weight_FMT_k_space_%s.txt", species[i].id);
                FILE* file_kernel = fopen(filename_fmt_kernel, "r");
                
                for (j = 0; j <N_kx; j++)
                {
                      for (l = 0; l < 3; l++)
                      {
                          fscanf(file_kernel, "%lf", &bdm);
                      }
                      for (l = 0; l < t_kernel_variables; l++)
                      {
                          fscanf(file_kernel, "%lf", &kernel_data[i][l][j][0]);
                      }
                }
                fclose(file_kernel);
            }
       }
       
        
       /*
       /// print the kernel data for the purpose of verification ... 
       for (i = 0; i< n_species; i++)
       {
            printf("The FMT kernel is given for the species %s : \n\n\n", species[i].id);
            for (int k = 0; k < N_kx; k++)
            {
               for (j=0; j < t_kernel_variables; j++)
               {    
                    printf(" %lf ", kernel_data[i][j][k][0]); // Real part
               }
               printf("\n");
            }
            
            printf("\n\n\n\n");
        
       }
        */
       /////////////////////////////////////////////////////////////////////////
       
       /////////////////////////////////////////////////////////////////////////
       /// input data is being fed for the rho values...
       FILE* file_rho = fopen("rho_mue_r_space.txt", "r");
       
       
       for (i = 0; i < N_x; i++ )
       {
            for (j = 0; j < n_species; j++)
            {
                  if (j == 0 )
                  {
                        fscanf(file_rho, "%lf", &r_data_fmt_1d[j][t_variables-1][i][0]);
                        fscanf(file_rho, "%lf", &bdm);
                        fscanf(file_rho, "%lf", &bdm);
                  }
                  else 
                  {
                        r_data_fmt_1d[j][t_variables-1][i][0] = r_data_fmt_1d[0][t_variables-1][i][0];
                  
                  }
                  fscanf(file_rho, "%lf", &r_data_fmt_1d[j][0][i][0]);
                  fscanf(file_rho, "%lf", &bdm);
            }
       }
       
       fclose(file_rho);
       
       
       /*
       for (i = 0; i < N_x; i++ )
       {
            for (j = 0; j < n_species; j++)
            {
                  printf("%lf   ", r_data_fmt_1d[j][0][i][0]);     
            }
            printf("\n");
       }
       */
       ///////////////////////////////////////////////////////////////////////// 
      
       /////////////////////////////////////////////////////////////////////////
       
       FILE* file_k_space = fopen("k_space.txt", "r");
       
       for (i = 0; i < N_kx; i++ )
       {
            for (j = 0; j < n_species; j++)
            {
                  if (j == 0 )
                  {
                        fscanf(file_k_space, "%lf", &k_data_fmt_1d[j][t_variables-1][i][0]);
                        fscanf(file_k_space, "%lf", &bdm);
                        fscanf(file_k_space, "%lf", &bdm);
                  }
                  else 
                  {
                        k_data_fmt_1d[j][t_variables-1][i][0] = k_data_fmt_1d[0][t_variables-1][i][0];
                  
                  }
            }
       }
       
       fclose(file_k_space);
       
       /*
       for (i = 0; i < N_x; i++ )
       {
            for (j = 0; j < n_species; j++)
            {
                  printf("%lf   ", k_data_fmt_1d[j][t_variables-1][i][0]);     
            }
            printf("\n");
       }
       */
       /////////////////////////////////////////////////////////////////////////
       
       ///###################################################################///
       
       
       
       
       
       
       
       
       
       ///###################################################################///
       
       
       /////////////////////////////////////////////////////////////////////////
       /// ... defining fftw_plans ..........................................///
       
       
       
       
       
       /////////////////////////////////////////////////////////////////////////
       ///... Picard's iteration sections is being implemented ................///
       
      
       
        fftw_plan* fftw_plan_forward_rho =  malloc(n_species * sizeof (fftw_plan));
        for (i = 0; i < n_species; i++)
        {
            fftw_plan_forward_rho[i] = fftw_plan_dft_1d(N_x, r_data_fmt_1d[i][0], k_data_fmt_1d[i][0], FFTW_FORWARD, FFTW_ESTIMATE);
        }
        
        
        fftw_plan** fftw_plan_backward_rho = malloc(n_species * sizeof (fftw_plan*));
        for (i = 0; i < n_species; i++)
        {
                fftw_plan_backward_rho[i] = malloc (6* sizeof(fftw_plan));
        }
        
        
        for (i = 0; i < n_species; i++)
        {
            for (j = 0; j < 6; j++)
            {
                
                fftw_plan_backward_rho[i][j] = fftw_plan_dft_1d(N_x, k_data_fmt_1d[i][j+13], r_data_fmt_1d[i][j+1], FFTW_BACKWARD, FFTW_ESTIMATE);
            }
          
        }
       
        fftw_plan** fftw_plan_forward_dphi = malloc(n_species * sizeof (fftw_plan*));
        for (i = 0; i < n_species; i++)
        {
            fftw_plan_forward_dphi[i] = malloc(6* sizeof(fftw_plan));
        }
        
         for (i = 0; i < n_species; i++)
        {
            for (j = 0; j < 6; j++)
            {
                fftw_plan_forward_dphi[i][j] = fftw_plan_dft_1d(N_x, r_data_fmt_1d[i][j+7], k_data_fmt_1d[i][j+7], FFTW_BACKWARD, FFTW_ESTIMATE);
            }
          
        }
        
        
        fftw_plan** fftw_plan_backward_dphi = malloc(n_species * sizeof (fftw_plan*));
        for (i = 0; i < n_species; i++)
        {
            fftw_plan_backward_dphi[i] = malloc(6* sizeof(fftw_plan));
        }
        
         for (i = 0; i < n_species; i++)
        {
            for (j = 0; j < 6; j++)
            {
                fftw_plan_backward_dphi[i][j] = fftw_plan_dft_1d(N_x, k_data_fmt_1d[i][j+13], r_data_fmt_1d[i][j+13], FFTW_BACKWARD, FFTW_ESTIMATE);
            }
          
        }
        
        
        
   
        
       
       
        calculate_rho_fmt_1d(k_data_fmt_1d, r_data_fmt_1d, kernel_data, fftw_plan_forward_rho, fftw_plan_backward_rho, fftw_plan_forward_dphi, fftw_plan_backward_dphi, N_kx, n_species, flag_fmt);
       
        
        printf("\n\n...here is the problem please consider it... \n\n");
       
        
        // Free allocated memory
        for ( i = 0; i < n_species; i++) {
            for ( j = 0; j < t_variables; j++) {
                free(k_data_fmt_1d[i][j]);
            }
            free(k_data_fmt_1d[i]);
        }
        free(k_data_fmt_1d);

        for ( i = 0; i < n_species; i++) {
            for ( j = 0; j < t_variables; j++) {
                free(r_data_fmt_1d[i][j]);
            }
            free(r_data_fmt_1d[i]);
        }
        free(r_data_fmt_1d);

        for ( i = 0; i < n_species; i++) {
            for ( j = 0; j < t_kernel_variables; j++) {
                free(kernel_data[i][j]);
            }
            free(kernel_data[i]);
        }
        free(kernel_data);
        
        
        
        
    
    }
    
    
    
    
  
    
    
    
    
    

    return 0;
}

