#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include <omp.h>
#include <fftw3.h>

const int E=20, UN=200, na=60;

long int tsim, ci, dpuc, dpuc1, ve, vel;




int  N, i, j, k, l, m, n;

double  rho_total, temperature, pressure, fugacity, bdm; 

double  current, adder, threshold;  



int main(void)
{

	piee=acos(-1.0);
	dis=piee/(float)na;
	tpiee=2.0*piee;
	EPSILON=0.000001;

        printf("\n\n\n\n\n\n****************************   Thanks for choosing me. I am going to minimize the given density functional in the NVT ensemble  !!!!! ******************************* \n\n\n");
        
        FILE* fu;
        fu = fopen("Data.txt", "r");
        
        
        
        
        
        
         // Parameter reading !!!!!!!!!!!!!!!
        
    
        
        int box type;
        
        fscanf(fu, "%lf", &bdm);                                     /// functional minimization would be done on this box type !!!!!!
        
        box_type=bdm;                                                  /// default box size is a Cartesian 3D box !!!!!! 
        
        double ble[3];                                               /// in Cartesian coordinate box size!!!!!!
        
        double shr[2];                                               /// in Curvy coordinate box size!!!!!!
        
        double dummy[3];                                             /// dummy size later supplied to the appropriate box dimensional variables e.g, size, shape etc. !!!!!! 
        
        dummy[0]=10.0;                                               /// default box is the 3D rectangular box and the default box size if 10!!!!!!!!!!!
        dummy[1]=10.0;
        dummy[2]=10.0;
        
        fscanf(fu, "%lf", &bdm);
        dummy[0]=bdm;
        fscanf(fu, "%lf", &bdm);
        dummy[1]=bdm;
        fscanf(fu, "%lf", &bdm);
        dummy[2]=bdm;
         
        
        if (box_type==0)
        {
            printf("No specification have ever been done as the box is chosen to be zero dimensional !!!!\n\n");
        }
        else if (box_type==1)
        {   
            ble[2]=dummy[2];
        }
        else if (box_type==2)
        {  
            ble[0]=dummy[0];
            ble[1]=dummy[1];
        }
        else if (box_type==3)
        {
            ble[0]=dummy[0];
            ble[1]=dummy[1];
            ble[2]=dummy[2];
        }
        else if (box_type==4)
        { 
            shr[0]=dummy[0];
            shr[1]=dummy[1];
        }
        else if (box_type==5)
        {
            shr[0]=dummy[0];
        }

    
        
        
        
        
        
        
        
        int war;                                                     /// type of particles  !!!!!!
        
        
        fscanf(fu, "%lf", &bdm);
        war=(int) ciel(bdm)
        
      
      
      
          
          
          
          
          
      
        int** inter_type=malloc(war*sizeof(int*));
        
        for (i=0; i<war; i++)
        {
            inter_type[i]=malloc(war*sizeof(int));
        }
        
        
        for (i=0; i<war; i++)
        
        {
            for(j=0; j<war; j++)
            {
                inter_type[i][j]=0;                                   /// setting default interaction which is just the hard-core repulsion
            }
        }
        
        
        fscanf(fu, "%lf", &bdm);
        for (i=0; i<war; i++)
        {
            for(j=i; j<war; j++)
            {
                fscanf(fu, "%lf", bdm);
                if(bdm==50.0) break;
                else   inter_type[i][j]=(int)bdm;
                
                inter_type[j][i]=inter_type[i][j];
                
            }   
        }
      
        
        for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
        
        
        
        
        
        
        
        
        
        
        double** ptc_size=malloc(war*sizeof(double*));
        
        for (i=0; i< war; i++)
        {
            ptc_size[i]=malloc(war*sizeof(double));
        }
        
        for (i=0; i<war; i++)
        {
            for (j=0; j<war; j++)
            {
                ptc_size[i][j]=1.0;                                   /// setting default particle size which is just the 1 unit !!!!!!
            }
        }
        
        fscanf(fu, "%lf", &bdm);
        for (i=0; i<war; i++)
        {
            for(j=i; j<war; j++)
            {
                fscanf(fu, "%lf", bdm);
                if(bdm==50.0) break;
                else   ptc_size[i][j]=bdm;
                ptc_size[j][i]=ptc_size[i][j];
                
            }   
        }
             
        
        for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
        
        
        
        
        
        
        
        
        
        
        
        double**  interaction_strength=malloc(war*sizeof(double*));
        
        for (i=0; i<war; i++)
        {
            interaction_strength[i]=malloc(war*sizeof(double));
        }
        
        for (i=0; i<war; i++)
        {
            for (j=0; j<war; j++)
            {
                interaction_strength[i][j]=-1.0;                      /// setting default interaction strength !!!!!!
            }
        }
        
        
        
        fscanf(fu, "%lf", &bdm);
        for (i=0; i<war; i++)
        {
            for(j=i; j<war; j++)
            {
                fscanf(fu, "%lf", bdm);
                if(bdm==50.0) break;
                else   interaction_strength[i][j]=bdm;
                interaction_strength[j][i]=interaction_strength[i][j];
                
            }   
        }
             
        
        for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
	
	
	
	
	
	
	double**  cutoff_range=malloc(war*sizeof(double*));
        
        for (i=0; i<war; i++)
        {
            cutoff_range[i]=malloc(war*sizeof(double));
        }
        
        for (i=0; i<war; i++)
        {
            for (j=0; j<war; j++)
            {
                cutoff_range[i][j]=3.0;                               /// setting default cutoff range !!!!!!
            }
        }
        
        
        
        fscanf(fu, "%lf", &bdm);
        for (i=0; i<war; i++)
        {
            for(j=i; j<war; j++)
            {
                fscanf(fu, "%lf", bdm);
                if(bdm==50.0) break;
                else   cutoff_range[i][j]=bdm;
                cutoff_range[j][i]=cutoff_range[i][j];
                
            }   
        }
             
        
        for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
        
        
        
        
        
      
        
        
        
        
        
        
                                                                     /// total density in the system !!!!!!
        
        fscanf(fu,"%lf", &bdm);
        rho_total = bdm;
        
    
        
        double* frac_ind=malloc(war*sizeof(double));                  /// fraction of each component in the system !!!!!!
        double* rho_ind=malloc(war*sizeof(double));                   /// density of the individual component in the system !!!!!!
        
        
        adder=1.0/war;
        current=0.0;
        
        for (i=0; i<war-1; i++)
        {
            frac_ind[i]=adder;
            current=current+adder;
        }
        
        frac_ind[war-1]=1.0-current; 
        
        
        fscanf(fu, "%lf", &bdm);
        for (i=0; i<war; i++)
        {
        
                fscanf(fu, "%lf", bdm);
                if(bdm==50.0) break;
                else   frac_ind[i]=bdm;
                
            
        }
        
        
        for (i=0; i<100; i++)
	{
		if (bdm==50.0) break;
		fscanf(bu,"%lf",&bdm);
	}
        
        
        for (i=0; i<war; i++)
        {
                rho_ind[i]= rho_total*frac_ind[i]
        }
        
        
      
        fscanf(fu, "%lf", &bdm);                                      /// setting temperature given to the system
        temp=bdm;
        
      
        
        fscanf(fu, "%lf", &bdm);                                      /// setting precision for the calculation
        dz=bdm;
        
        
        
        
        
        
        
        
        
        
        
        /// Thanks for the parameter specification, now the real calculation will start !!!!!!
        
        
        
        
        
        
        
        
        if (box_type==1)
        {
                                                                       /// Allocate memory for input and output arrays
              
                                                                       /// For the simplicity of the problems the variables are defined for each component in the full set 
               
                        int N_z;                                       /// number of points in the real space
                        
                        
                        N_z = (int) (ble[2]/dz) +1;
                      
                        
                        fftw_complex** rho_k_3 = fftw_malloc(war * sizeof(fftw_complex*));
                        fftw_complex** rho_z = fftw_malloc(war*sizeof(fftw_complex*));
                        fftw_complex*** new_rho_z = fftw_malloc(6 * sizeof(fftw_complex**));
                        
                        
                        for (i=0; i<war; i++)
                        {
                              rho_k_3[i] = fftw_malloc(N_z * sizeof(fftw_complex));
                              rho_z[i] = fftw_malloc(N_z * sizeof(fftw_complex));
                        }
                        
                        for (i=0; i<6; i++)
                        {
                              new_rho_z[i]=fftw_malloc(war*sizeof(fftw_complex*))
                              for (j=0; j<war; j++)
                              {
                                  new_rho_z[i][j]=fftw_malloc(N_z*sizeof(fftw_complex));
                              }
                        }
                        
                        
                        fftw_plan* forward_plan = malloc(war* sizeof(fftw_plan));
                        
                        
                        fftw_plan** backward_plan = malloc(6* sizeof(fftw_plan*));
                        
                        for(i=0; i<6; i++)              
                        (
                                backward_plan[i]= malloc(war*sizeof(fftw_plan));
                        }
                        
                        
                        

                        
                        double*** FMT_kernel_3=malloc(6*sizeof(double**));
                        fftw_complex*** FMT_rhok_kernel_3=malloc(6*sizeof(fftw_complex**));
                        
                        for (i=0; i<6 ; i++)
                        { 
                              FMT_kernel_3[i]=malloc(war*sizeof(double*));
                              FMT_rhok_kernel_3[i]=malloc(war*sizeof(fftw_complex*));
                              
                              for(j=0; j<war; j++)
                              {
                                    FMT_kernel_3[i][j]=malloc(N_z*sizeof(double));
                                    FMT_rhok_kernel_3[i][j]=malloc(N_z*sizeof(fftw_complex));
                              }
                              
                        
                        }
                        
                        
                        
                        
                        
                        
                        
                        for (i=0; i<war; i++)
                        {
                                
                                forward_plan[i] = fftw_plan_dft_1d(N_z, rho_z[i], rho_k_3[i], FFTW_FORWARD, FFTW_ESTIMATE);
                                
                 
                        }
                        

                        
                        for (i=0; i<6; i++)
                        {
                            for(j=0; j<war; j++)
                            {
                                backward_plan[i][j]= fftw_plan_dft_1d(N_z, FMT_rhok_kernel_3[i][j], new_rho_z[i][j], FFTW_BACKWARD, FFTW_ESTIMATE); 
                            }
                        }

                                                                           // to create forward and backward FFT plans
                        forward_plan = fftw_plan_dft_1d(N_z, rho_z, rho_k_z, FFTW_FORWARD, FFTW_ESTIMATE);
                        backward_plan = fftw_plan_dft_1d(N_z, rho_k_3, new_rho_z, FFTW_BACKWARD, FFTW_ESTIMATE);

                                                                           // to initialize density profile rho (initial guess)
                        for (i = 0; i < war; i++) 
                        {
                            for(j=0; j<N_z; j++)
                            {
                                rho_z[i][j][0] = rho_ind[i]; // Real part (density guess)
                                rho_z[i][j][1] = 0.0;                      // imaginary part (zero for real input)
                            }
                        }

                     
                        
                        

                        // Picard iteration loop
                        
                        
                        double* kspace=malloc(N_z*sizeof(double));
                        double* rspace=malloc(N_z*sizeof(double));
                        
                        for(i=0; i<N_z; i++)
                        {
                              
                              rspace=dz*i;
                        }
                        
                      
                        
                        // Compute the wavevector frequencies
                        
                        for (int i = 0; i < N_z; i++) 
                        {
                            if (i <= N_z / 2) 
                            {
                                kspace[i] = i * (tpiee / ((N_z-1) * d));  // Positive frequencies
                            } 
                            else 
                            {
                                kspace[i] = (i - N_z) * (tpiee / ((N_z-1) * d));  // Negative frequencies
                            }
                        }
                        
                        
                        
                        
                        
                                                                      /// kernel in k-space
                                                                      
                                                                      
                        for (i=0; i<N_z; i++)
                        {
                            if (kspace<EPSILON) 
                            {
                                  for (j=0; j<war ; j++)
                                  {
                                      FMT_kernel_3[0][j][i]=1;
                                      FMT_kernel_3[1][j][i]=ptc_size[j][j]*0.5;
                                      FMT_kernel_3[2][j][i]=piee*ptc_size[j][j]*ptc_size[j][j]*0.25;
                                      FMT_kernel_3[3][j][i]=piee*0.125*ptc_size[j][j]*ptc_size[j][j]*ptc_size[j][j];
                                      FMT_kernel_3[4][j][i]=0;
                                      FMT_kernel_3[5][j][i]=0;
                                  }
                                  
                            }
                            else    
                            {
                                  for (j=0; j<war; j++)
                                  {
                                      FMT_kernel_3[0][j][i] = sin(kspace[i]*piee*ptc_size[j][j])/(kspace[i]*ptc_size[j][j]*piee) ; /// alpha==0
                                      
                                      FMT_kernel_3[1][j][i] = sin(kspace[i]*piee*ptc_size[j][j])/(2.0*k_space[i]*piee)  ;/// alpha==1
                                    
                                      FMT_kernel_3[2][j][i] = ptc_size[j][j]*sin(kspace[i]*piee*ptc_size[j][j])/kspace[i]  ;/// alpha==2
                                      
                                      FMT_kernel_3[3][j][i] = sin(kspace[i]*piee*ptc_size[j][j])/(2.0*pow(kspace[i], 3.0)*piee*piee)   - ptc_size[j][j] * cos(kspace[i]*piee*ptc_size[j][j])/ (2.0* kspace[i]*kspace[i]*piee);/// alpha==3
                                      
                                      FMT_kernel_3[4][j][i] = (kspace[i]*piee*ptc_size[j][j]*cos(kspace[i]*piee*ptc_size[j][j])-sin(kspace[i]*piee*ptc_size[j][j]))/(2.0*ptc_size[j][j] *piee*piee*kspace[i]*kspace[i]);  /// alpha==4
                                    
                                      FMT_kernel_3[5][j][i]=  (kspace[i]*piee*ptc_size[j][j]*cos(kspace[i]*piee*ptc_size[j][j])-sin(kspace[i]*piee*ptc_size[j][j]))/(kspace[i]*kspace[i]*piee);
                                      
                                  }
                        
                            }
                        }
                        
                        
                        
                        
                        
                        
                        
                        
                        
                            
                            
                        fftw_complex** d_phi=malloc(6*sizeof(fftw_complex*));
                        fftw_complex** d_phi_k_3=malloc(6*sizeof(fftw_complex*));
                        fftw_complex*** new_d_phi=malloc(6*sizeof(fftw_complex**));
                        
                        fftw_complex*** kernel_d_phi_k3=malloc(6*sizeof(fftw_complex**));
                        
                       
                        for(i=0; i<6; i++)
                        { 
                            d_phi[i]=malloc(N_z*sizeof(fftw_complex));  
                            d_phi_k_3[i]=malloc(N_z*sizeof(fftw_complex));  
                            new_d_phi[i]=malloc(war*sizeof(fftw_complex*));
                            kernel_d_phi_k3[i]=malloc(war*sizeof(fftw_complex*));
                            for (j=0; j<war; j++)
                            {
                                new_d_phi[i][j]=malloc(N_z*sizeof(fftw_complex))
                                kernel_d_phi_k3[i][j]=malloc(N_z*sizeof(fftw_complex))
                            }
                        }
                        
                        for (i=0; i<6; i++)
                        {
                            for (j=0; j<N_z; j++)
                            {
                                d_phi[i][j][0]=0.0;
                                d_phi[i][j][1]=0.0;
                                d_phi_k_3[i][j][0]=0.0;
                                d_phi_k_3[i][j][1]=0.0; 
                                for (k=0; k<war; k++)
                                {
                                    new_d_phi[i][k][j][0]=0.0;
                                    new_d_phi[i][k][j][1]=0.0;
                                    kernel_d_phi_k3[i][k][j][0]=0.0;
                                    kernel_d_phi_k3[i][k][j][1]=0.0;
                                }
                            }                        
                        }
                        
                        
                        
                        
                        
                        
                        fftw_plan* forward_plan_phi = malloc(6* sizeof(fftw_plan));
                        
                        for (i=0; i<6; i++)
                        {
                                
                                forward_plan_phi[i] = fftw_plan_dft_1d(N_z, d_phi[i], d_phi_k_3[i], FFTW_FORWARD, FFTW_ESTIMATE);
                                
                 
                        }
                        
                        
                        
                        
                        
                        
                        
                        fftw_plan** backward_plan_phi = malloc(6* sizeof(fftw_plan*));
                        
                        for(i=0; i<6; i++)              
                        (
                                backward_plan_phi[i]= malloc(war*sizeof(fftw_plan));
                        }
                        
                        for (i=0; i <6; i++)
                        {
                            for (j=0; j<war; j++)
                            {
                                backward_plan_phi[i][j]=fftw_plan_dft_1d(N_z, kernel_d_phi_k3[i][j], new_d_phi[i][j], FFTW_BACKWARD, FFTW_ESTIMATE); 
                            }
                            
                        }
                        
                        
                        fftw_complex* total_n_aplha=malloc(6*sizeof(fftw_complex));
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        //// Dispersion variables !!!!!!!!!!
                        
                        fftw_complex*** uk_3=malloc(war*sizeof(fftw_complex**));
                        
                        for (i=0; i <war; i++)
                        {
                            uk_3=malloc(war*sizeof(fftw_complex*));
                            for (j=0; j<war. j++)
                            {
                                uk_3[i][j]=malloc(N_z*sizeof(fftw_complex))
                            }
                            
                        }
                        
                        
                        for (i=0; i<war; i++)
                        {
                            for (j=0; j<war; j++)
                            {
                                for(k=0; k<N_z; k++)
                                {
                                    uk_3[i][j][k][0]=0.0;
                                    uk_3[i][j][k][1]=0.0;
                                }
                            }
                        }
                        
                        
                        
                        
                        for (i=0; i<war; i<++)
                        {
                            for (j=0; j<war; j++)
                            {
                                
                                if(inter_type==1)
                                {
                                    for (k=0; k<N_z; k++)
                                    {
                                        uk_3[i][j][k][0]=interaction_strength[i][j]*pow(piee, 0.5)*ptc_size[i][j]*exp(-k_space[k]*k_space[k]*ptc_size[i][j]*ptc_size[i][j]/4.0);       
                                    }
                                else
                                {
                                    for (k=0; k<N_z; k++)
                                    {
                                        uk_3[i][j][k][0]=0.0;
                                    }
                                }
                            }
                        }
                        
                        
                        
                        
                        
                        
                        fftw_complex*** uk_rhok_3=malloc(war*sizeof(fftw_complex**));
                        
                        for(i=0; i<war; i++)
                        {
                            uk_rhok_3[i]=malloc(war*sizeof(fftw_complex*));
                            for(j=0; j<war; j++)
                            {
                                uk_rhok_3[i][j]=malloc(N_z*sizeof(fftw_complex));
                            }
                        }
                        
                        
                        for (i=0; i<war; i++)
                        {
                            for (j=0; j<war; j++)
                            {
                                for(k=0; k<N_z; k++)
                                {
                                    uk_rhok_3[i][j][k][0]=0.0;
                                    uk_rhok_3[i][j][k][1]=0.0;
                                }
                            }
                        }
                        
                        
                        
                        
                        
                        fftw_complex***  ddispersion=malloc(war*sizeof(fftw_complex**));
                        
                        for(i=0; i<war; i++)
                        {
                            ddispersion[i]=malloc(war*sizeof(fftw_complex*));
                            for(j=0; j<war; j++)
                            {
                                ddispersion[i][j]=malloc(N_z*sizeof(fftw_complex));
                            }
                        }

                        
                        for (i=0; i<war; i++)
                        {
                            for (j=0; j<war; j++)
                            {
                                for (k=0; k<N_z; k++)
                                {
                                    ddispersion[i][j][k][0]=0.0;
                                    ddispersion[i][j][k][1]=0.0;
                                }
                            }
                        }
                        
                        
                        
                        
                        fftw_plan** backward_plan_disp = malloc(war* sizeof(fftw_plan*));
                        
                        for(i=0; i<war; i++)              
                        (
                                backward_plan_disp[i]= malloc(war*sizeof(fftw_plan));
                        }
                        
                        for (i=0; i <war; i++)
                        {
                            for (j=0; j<war; j++)
                            {
                                backward_plan_disp[i][j]=fftw_plan_dft_1d(N_z, uk_rhok_3[i][j], ddispersion[i][j], FFTW_BACKWARD, FFTW_ESTIMATE); 
                            }
                            
                        }
                        
                        
                        fftw_complex* total_n_aplha=malloc(6*sizeof(fftw_complex));
                        
                        
                        
                        int max_iterations = 1000;
                        double tolerance = 1e-6;
                        double error = 1.0;
                        int iter = 0;
                        
                        
                        while (iter < max_iterations && error > tolerance) 
                        {
                                    ///              FMT     parameters  !!!!!!!!!!!!
                                    
                                    
                                    
                                    
                                    for ( i = 0; i < war; i++) 
                                    {    
                                        fftw_execute(forward_plan[i]);  // FFT on each component
                                    }

                                    
                                    
                                    
                                    
                                    for (i=0; i<4; i++)
                                    {
                                        for (j=0; j<war; j++)
                                        
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   FMT_rhok_kernel_3[i][j][k][0]=rho_k_3[j][k][0]*FMT_kernel_3[i][j][k];
                                                   FMT_rhok_kernel_3[i][j][k][1]=rho_k_3[j][k][1]*FMT_kernel_3[i][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                    }
                                    
                                        for (j=0; j<war; j++)
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   FMT_rhok_kernel_3[4][j][k][0]=-rho_k_3[j][k][1]*FMT_kernel_3[4][j][k];
                                                   FMT_rhok_kernel_3[4][j][k][1]=rho_k_3[j][k][0]*FMT_kernel_3[4][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                        
                                        
                                        for (j=0; j<war; j++)
                                        
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   FMT_rhok_kernel_3[5][j][k][0]=-rho_k_3[j][k][1]*FMT_kernel_3[5][j][k];
                                                   FMT_rhok_kernel_3[5][j][k][1]=rho_k_3[j][k][0]*FMT_kernel_3[5][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                    
                                
                                    
                                    for (i=0; i<6; i++)
                                    {
                                        for(j=0; j<war; j++)
                                        {
                                            fftw_exicute(backward_plan[i][j]) ;    ////    = fftw_plan_dft_1d(N_z, FMT_rhok_kernel_3[i][j], new_rho_z[i][j], FFTW_BACKWARD, FFTW_ESTIMATE); 
                                        }
                                    }

                                    
                                    
                                    
                                    ///   FMT calculation !!!!!!
                                    
                                    
                                                                    /// calculating dphi/dalpha   !!!!!!
                                    
                                    
                                    for (j=0; j <N_z; j++)
                                    {
                                          for (i=0; i<6; i++)
                                          {
                                              total_n_alpha[i][0]=0.0;
                                              total_n_alpha[i][1]=0.0;
                                               
                                              for(k=0; k<war; k++)
                                              {
                                                  total_n_alpha[i][0]=total_n_alpha[i][0]+new_rho_z[i][k][j][0];
                                                  total_n_alpha[i][1]=total_n_alpha[i][0]+new_rho_z[i][k][j][1];
                                              }   
                                              
                                              d_phi[i][j][1]=0.0;
                                        
                                          } 
                                          
                                          d_phi[0][j][0]= -log(1.0-total_n_alpha[3][0]);
                                          
                                          d_phi[1][j][0]= total_n_alpha[2][0]/(1.0-total_n_alpha[3][0]);
                                          
                                          d_phi[2][j][0]= total_n_alpha[1][0]/(1.0-total_n_alpha[3][0]) + (1.0/(8.0*piee))*(total_n_alpha[2][0]*total_n_alpha[2][0]- total_n_alpha[4][0]*total_n_alpha[4][0]) /((1.0-total_n_alpha[3][0]) * (1.0-total_n_alpha[3][0]));
                                          
                                          d_phi[3][j][0]= total_n_alpha[0][0]/(1.0-total_n_alpha[3][0])  +  (total_n_alpha[1][0] * total_n_alpha[2][0]-total_n_alpha[4][0]*total_n_alpha[5][0])/((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0])) ;
                                          
                                          d_phi[3][j][0]= d_phi[3][j][0]+(1.0/(12.0*piee))*(total_n_alpha[2][0]*total_n_alpha[2][0]*total_n_alpha[2][0]-3.0*total_n_alpha[5][0]*total_n_alpha[5][0])/ ((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0]));
                                          
                                          d_phi[4][j][0]= total_n_alpha[5][0]/(1.0-total_n_alpha[3][0]);
                                          
                                          d_phi[5][j][0]= total_n_alpha[4][0]/(1.0-total_n_alpha[3][0]) + 0.25  * (total_n_alpha[2][0]*total_n_alpha[5][0])/((1.0-total_n_alpha[3][0])*(1.0-total_n_alpha[3][0]));
                                          
                                          
                                          
                                    }
                                   
                                   
                                   
                                   
                                    
                                    
                                    
                                    
                                    
                                    for (i=0; i<6; i++)
                                    {
                                        fftw_exicute(forward_plan_phi[i]);
                                    }                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                     for (i=0; i<4; i++)
                                    {
                                        for (j=0; j<war; j++)
                                        
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   kernel_d_phi_k3[i][j][k][0]=d_phi_k_3[i][k][0]*FMT_kernel_3[i][j][k];
                                                   kernel_d_phi_k3[i][j][k][1]=d_phi_k_3[i][k][1]*FMT_kernel_3[i][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                    }
                                    
                                        for (j=0; j<war; j++)
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   kernel_d_phi_k3[4][j][k][0]=-d_phi_k_3[4][k][1]*FMT_kernel_3[4][j][k];
                                                   kernel_d_phi_k3[4][j][k][1]=d_phi_k_3[4][k][0]*FMT_kernel_3[4][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                        
                                        
                                        for (j=0; j<war; j++)
                                        
                                        {
                                            for (k=0 ; k<N_z; k++)
                                            {
                                                   kernel_d_phi_k3[5][j][k][0]=-d_phi_k_3[5][k][1]*FMT_kernel_3[5][j][k];
                                                   kernel_d_phi_k3[5][j][k][1]=d_phi_k_3[5][k][0]*FMT_kernel_3[5][j][k];
                                            
                                            }
                                            
                                        
                                        }
                                    
                                    
                                    
              
                                    
                                    for (i=0; i<6; i ++)
                                    {
                                        for (j=0; j<war; j++)
                                        {
                                            fftw_exicute(backward_plan_phi[i][j]);
                                        
                                        }
                                    
                                    }
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    /// Dispersion Implementation  .....  !!!!!!
                                    
                                    
                                    
                                    for (i=0; i<war; i++)
                                    {
                                        for (j=i; j<war; j++)
                                        {
                                            for (k=0; k<N_z; k++)
                                            {
                                                uk_rhok_3[i][j][k][0]=uk_3[i][j][k][0]*rhok_3[i][j][k][0];
                                                uk_rhok_3[i][j][k][1]=0.0;
                                            }
                                        }
                                    }
                                    
                                  
                                    for (i=0; i<war; i++)
                                    {
                                        for (j=i; j<war; j++)
                                        {
                                            fftw_exicute(backward_plan_dispersion[i][j]);   
                                        }
                                    }
                                    
                                    
                                    
                                    //// Mean field treatment !!!!!!
                                    
                                    
                                    //// Dispersion part which is being done by using the mean field treatment !!!!!!
                                    
                                    
                                    
                                    
                                    

                                    // Step 2: Apply equations in Fourier space (for example, solving a convolution)
                                    // Modify `rho_k` based on FMT or DFT theory (this is model-dependent)
                                    // For example, you might update rho_k by multiplying with some kernel:
                                    for (int i = 0; i < N; i++) {
                                        double kernel = exp(-pow(i, 2));  // Example kernel function
                                        rho_k[i][0] *= kernel;  // Real part update
                                        rho_k[i][1] *= kernel;  // Imaginary part update
                                    }

                                    // Step 3: Perform inverse FFT to return to real space
                                    fftw_execute(backward_plan);

                                    // Normalize the inverse FFT result
                                    for (int i = 0; i < N; i++) {
                                        new_rho[i][0] /= N;  // Real part (density update)
                                        new_rho[i][1] /= N;  // Imaginary part
                                    }

                                    // Step 4: Update density using Picard iteration (linear mixing)
                                    error = 0.0;
                                    for (int i = 0; i < N; i++) {
                                        double old_rho_real = rho[i][0];
                                        double old_rho_imag = rho[i][1];
                                        
                                        // Linear mixing: new_rho = alpha * new_rho + (1 - alpha) * old_rho
                                        double alpha = 0.5; // Mixing parameter
                                        rho[i][0] = alpha * new_rho[i][0] + (1 - alpha) * old_rho_real;
                                        rho[i][1] = alpha * new_rho[i][1] + (1 - alpha) * old_rho_imag;
                                        
                                        // Calculate error (e.g., L2 norm of difference between old and new densities)
                                        error += pow(rho[i][0] - old_rho_real, 2) + pow(rho[i][1] - old_rho_imag, 2);
                                    }

                                    error = sqrt(error);  // Root mean square error

                                    // Print iteration info
                                    printf("Iteration %d, error: %f\n", iter, error);

                                    iter++;
                        }

                        
              
        
        
        
        }
        
        
        
    
        






        // Allocate memory for input and output arrays
        
        fftw_plan forward_plan, backward_plan;
        
        
        fftw_complex* rho_k = fftw_malloc(N * sizeof(fftw_complex));
        fftw_complex* rho = fftw_malloc(N * sizeof(fftw_complex));
        fftw_complex* new_rho = fftw_malloc(N * sizeof(fftw_complex));
        
        
        
        

        // Create forward and backward FFT plans
        forward_plan = fftw_plan_dft_1d(N, rho, rho_k, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_1d(N, rho_k, new_rho, FFTW_BACKWARD, FFTW_ESTIMATE);

        // Initialize density profile rho (initial guess)
        for (int i = 0; i < N; i++) {
            rho[i][0] = sin(2 * M_PI * i / N); // Real part (density guess)
            rho[i][1] = 0.0;                    // Imaginary part (zero for real input)
        }

        int max_iterations = 1000;
        double tolerance = 1e-6;
        double error = 1.0;
        int iter = 0;

        // Picard iteration loop
        while (iter < max_iterations && error > tolerance) {
            // Step 1: Perform forward FFT to get rho in Fourier space (rho_k)
            fftw_execute(forward_plan);

            // Step 2: Apply equations in Fourier space (for example, solving a convolution)
            // Modify `rho_k` based on FMT or DFT theory (this is model-dependent)
            // For example, you might update rho_k by multiplying with some kernel:
            for (int i = 0; i < N; i++) {
                double kernel = exp(-pow(i, 2));  // Example kernel function
                rho_k[i][0] *= kernel;  // Real part update
                rho_k[i][1] *= kernel;  // Imaginary part update
            }

            // Step 3: Perform inverse FFT to return to real space
            fftw_execute(backward_plan);

            // Normalize the inverse FFT result
            for (int i = 0; i < N; i++) {
                new_rho[i][0] /= N;  // Real part (density update)
                new_rho[i][1] /= N;  // Imaginary part
            }

            // Step 4: Update density using Picard iteration (linear mixing)
            error = 0.0;
            for (int i = 0; i < N; i++) {
                double old_rho_real = rho[i][0];
                double old_rho_imag = rho[i][1];
                
                // Linear mixing: new_rho = alpha * new_rho + (1 - alpha) * old_rho
                double alpha = 0.5; // Mixing parameter
                rho[i][0] = alpha * new_rho[i][0] + (1 - alpha) * old_rho_real;
                rho[i][1] = alpha * new_rho[i][1] + (1 - alpha) * old_rho_imag;
                
                // Calculate error (e.g., L2 norm of difference between old and new densities)
                error += pow(rho[i][0] - old_rho_real, 2) + pow(rho[i][1] - old_rho_imag, 2);
            }

            error = sqrt(error);  // Root mean square error

            // Print iteration info
            printf("Iteration %d, error: %f\n", iter, error);

            iter++;
        }

        // Cleanup
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_cleanup();

        
        
        
        
        
        
        
        
        
        
        
        
        return 0;
        
        
        
        
        
        
        
        
        
        
        
        
        
}

        
        
        
        
        
        
        
        
        
        
        
        
        
       
}

        
        
        
        
        
        
        
      
        
        
        
        
        
        
       
        
         
        
        
          
        

      
	// Scaling the parameteres
	
	
	
	
	
	
	
}
