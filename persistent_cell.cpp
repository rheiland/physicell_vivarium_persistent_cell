
// https://docs.google.com/document/d/1rxG1U6g_l0XT-2-CtrJEFT8HypmM9ygQYfiNKjdYAZI/edit?tab=t.0#heading=h.ogpvfgg7xadq

// g++ persistent_cell.cpp -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11 -o persistent_cell

#include <iostream>
#include <random>
#include <omp.h>
#include <cstdio>

thread_local std::mt19937_64 physicell_PRNG_generator; 
thread_local bool local_pnrg_setup_done = false; 

// unsigned int physicell_random_seed = 0; 
std::vector<unsigned int> physicell_random_seeds; 

// std::vector<double> cell_xpos; 
// std::vector<double> cell_ypos; 

// double dt_mech = 0.1;
double max_time = 10000;   // mins
// double dt_mech = 0.5;    // default is 0.1
double dt_mech = 1.0;    // default is 0.1
std::vector<double> motility_vector {1.0, 0.0, 0.0};
std::vector<double> position = {50.,50.,0.};
std::vector<double> velocity = {0.1,0.,0.};
std::vector<double> previous_velocity = {0.1,0.,0.};

double UniformRandom( void )
{
	thread_local std::uniform_real_distribution<double> distribution(0.0,1.0);
	if( local_pnrg_setup_done == false )
	{
		// get my thread number 
		int i = omp_get_thread_num(); 
        // std::cout << "UniformRandom():  i (thread #) = " << i << std::endl;  // = 0
    	// physicell_PRNG_generator.seed( physicell_random_seeds[i] ); 
        int myseed = 42;
        // std::cout << "- 0)UniformRandom():  myseed= " << myseed << std::endl;
    	physicell_PRNG_generator.seed( myseed ); 
		local_pnrg_setup_done = true; 
/*
		#pragma omp critical 
		{
		std::cout << "thread: " << i 
		<< " seed: " << physicell_random_seeds[i]  << std::endl; 
		std::cout << "\t first call: " << distribution(physicell_PRNG_generator) << std::endl; 
		}
*/
	}
    double retval = distribution(physicell_PRNG_generator);
    // std::cout << "--1) UniformRandom():  retval= " << retval << std::endl;  // = 0
    // return distribution(physicell_PRNG_generator);
    return retval;

	// helpful info: https://stackoverflow.com/a/29710970
/*

	static std::uniform_real_distribution<double> distribution(0.0,1.0); 
	double out;
	out = distribution(physicell_PRNG_generator);
	return out; 
*/	
}

// In P*_utilities.cpp
std::vector<double> UniformOnUnitCircle( void )
{
	std::vector<double> output = {0,0,0}; 

	static long double two_pi = 6.283185307179586476925286766559;  
	                       
	long double theta = UniformRandom(); //  BioFVM::uniform_random();
	theta *= two_pi; // Choose theta uniformly distributed on [0, 2*pi).

	output[0] = cos(theta); 
	output[1] = sin(theta); // (cos(t) , sin(t) , 0 )

	return output; 
}

// BioFVM_vector.cpp
void normalize( std::vector<double>* v )
{
 double norm = 1e-32; 

 for( unsigned int i=0; i < v->size(); i++ )
 { norm += ( (*v)[i] * (*v)[i] ); }
 norm = sqrt( norm ); 

 for( unsigned int i=0; i < v->size(); i++ )
 { (*v)[i] /=  norm ; }

 // If the norm is small, normalizing doens't make sense. 
 // Just set the entire vector to zero. 
 static bool I_warned_you = false; 
 if( norm <= 1e-16 )
 { 
  if( I_warned_you == false )
  {
   std::cout << "Warning and FYI: Very small vectors are normalized to 0 vector" << std::endl << std::endl; 
   I_warned_you = true; 
  }

  for( unsigned int i=0; i < v->size(); i++ )
  { (*v)[i] = 0.0; }
 }

 return; 
}

// BioFVM_vector.cpp
void axpy( std::vector<double>* y, const double& a , const std::vector<double>& x )
{
 for( unsigned int i=0; i < (*y).size() ; i++ )
 {
  (*y)[i] += a * x[i] ; 
 }
 return ; 
}

// in core/P*_cell.cpp
//void Cell::update_motility_vector( double dt_ )
void update_motility_vector( double dt_ )
{
    std::vector<double> randvec(3,0.0);
    static double migration_bias = 0.5;
    // static double migration_bias = 1.0;
    static double migration_speed = 0.1;
    static double persistence_time = 0.0;

    // std::cout << "\n---- update_motility_vector():\n";

    // if( phenotype.motility.is_motile == false )
    // {
    //         phenotype.motility.motility_vector.assign( 3, 0.0 );
    //         return;
    // }
    // if( UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
    // std::cout << "---- update_motility_vector():  dt=" << dt_ << ", motility_vector=" << motility_vector[0] <<", "<< motility_vector[1] << std::endl;

    if( UniformRandom() < dt_ / persistence_time || persistence_time < dt_ )
    {
            // if( phenotype.motility.restrict_to_2D == true )
            { randvec = UniformOnUnitCircle(); 
            //   std::cout << "randvec= " << randvec[0] << ", " << randvec[1] << std::endl; 
            }
            // else
            // { randvec = UniformOnUnitSphere(); }
            // // if the update_bias_vector function is set, use it
            // if( functions.update_migration_bias )
            // {
            //         functions.update_migration_bias( this,phenotype,dt_ );
            // }
        
        // recall in custom.cpp, we did:
        // std::vector<double> ctype1_direction {1.0, 0.0, 0.0};
        // cell_defaults.phenotype.motility.migration_bias_direction = ctype1_direction;	


            // phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motility = bias_vector
            motility_vector = {1.,0.,0.};

            // std::cout << "motility_vector(2)= "<< motility_vector[0] << ", " << motility_vector[1] <<", "<< motility_vector[2]<< std::endl;

            // std::vector<double> motility_vector {1.0, 0.0, 0.0};   // rwh: do we really reset?
            // phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector
            // motility_vector *= migration_bias; // motility = bias*bias_vector
            motility_vector[0] *= migration_bias; // motility = bias*bias_vector
            motility_vector[1] *= migration_bias; // motility = bias*bias_vector
            // std::cout << "motility_vector (3)= "<< motility_vector[0] << ", " << motility_vector[1] << std::endl;
            // double one_minus_bias = 1.0 - phenotype.motility.migration_bias;
            double one_minus_bias = 1.0 - migration_bias;
            // std::cout << "one_minus_bias = "<< one_minus_bias << std::endl;
            // axpy( &(phenotype.motility.motility_vector), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
            axpy( &(motility_vector), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
            // std::cout << "motility_vector (4)= "<< motility_vector[0] << ", " << motility_vector[1] << std::endl;
            // for (size_t i = 0; i < motility_vector.size(); ++i) 
            // {
            //     randvec[i] += one_minus_bias * motility_vector[i];
            // }
            // normalize( &(phenotype.motility.motility_vector) );
            normalize( &(motility_vector) );
            // std::cout << "motility_vector (norm'd)= "<< motility_vector[0] << ", " << motility_vector[1] << std::endl;
            // phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
            // motility_vector *= migration_speed;
            motility_vector[0] *= migration_speed;
            motility_vector[1] *= migration_speed;
            // std::cout << "motility_vector (5)= "<< motility_vector[0] << ", " << motility_vector[1] << std::endl;

            // cell_xpos.push_back(motility_vector[0]);
            // cell_ypos.push_back(motility_vector[1]);

            // std::cout << "motility_vector (*speed)= "<< motility_vector[0] << ", " << motility_vector[1] << std::endl;
    }
    return;
}

// void Cell::update_position( double dt )
void update_position( double dt )
{
	// BioFVM Basic_Agent::update_position(dt) returns without doing anything. 
	// So we remove this to avoid any future surprises. 
	// 
	// Basic_Agent::update_position(dt);
		
	// use Adams-Bashforth 
	static double d1; 
	static double d2; 
	static bool constants_defined = false; 
	if( constants_defined == false )
	{
		d1 = dt; 
		d1 *= 1.5; 
		d2 = dt; 
		d2 *= -0.5; 
		constants_defined = true; 
	}
	
	// new AUgust 2017
	// if( default_microenvironment_options.simulate_2D == true )
	{ velocity[2] = 0.0; }
	
	// std::vector<double> old_position(position);   // rwh - do a PR to comment out this!
    // std::cout << "--- pre-update position:  d1= "<<d1<< ", d2=" << d2 << std::endl;
    // std::cout << "--- velocity:  "<<velocity[0]<<", "<<velocity[1]  << std::endl;
    // std::cout << "--- previous_velocity:  "<<previous_velocity[0]<<", "<<previous_velocity[1]  << std::endl;
	axpy( &position , d1 , velocity );           //  position += d1 * velocity
	axpy( &position , d2 , previous_velocity );  //  position += d1 * previous_velocity
	// overwrite previous_velocity for future use 
	// if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
		// std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;
	
	previous_velocity = velocity; 
	
	velocity[0]=0; velocity[1]=0; velocity[2]=0;
	// if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
	// {
	// 	updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	// }
	// else
	// {
	// 	updated_current_mechanics_voxel_index=-1;
		
	// 	is_out_of_domain = true; 
	// 	is_active = false; 
	// 	is_movable = false; 
	// }
	return; 
}


int main()
{
    omp_set_num_threads(1);
    // std::cout <<"cell pos= "<<position[0]<<","<<position[1]<<std::endl;

    std::cout<<"time,id,com_1,com_2,area,surface\n";
    // std::cout<<"time,id,com_1,com_2\n";
    // printf("time,id,com_1,com_2\n");

    int max_steps = int(max_time / dt_mech);
    int idx_output = int(100/dt_mech);
    // printf("dt_mech = %f\n",dt_mech);
    // printf("max_time = %f\n",max_time);
    // printf("max_steps = %d\n",max_steps);
    // printf("idx_output = %d\n",idx_output);

    int id_run = 0;
    for (int run_num=0; run_num<100; run_num++)
    // for (int id_run=100; id_run<=1000; id_run+=100)   
    {
        // time,id,com_1,com_2,area,surface
        // 100,0,56.3445,50.1373,222.34,52.86
        // 200,0,62.7075,50.2596,222.34,52.86

        // reset
        position[0] = 50.0;
        position[1] = 50.0;

        // Not sure where, but somewhere in the full PhysiCell sim, this is called N times
        // before we get into the actual motility. Inserting these additional UniformRandom
        // call would get us bitwise reprod (assuming same seed and 1 core).
        // if (idx==0)
        {
            // std::cout << " >>> extra UniformRandom calls\n";
            UniformRandom();
            UniformRandom();
            UniformRandom();
            UniformRandom();
        }

        // for (int idx=0; idx<100; idx++)   // (100 / dt_mech) = 1000
        // std::cout <<"0,"<< run_num<<","<< position[0]<<","<<position[1] << std::endl;
        // for (int idx=1; idx<=10; idx++)   // equiv of 1 sec: (1/ 0.1) = 10
        for (int idx=1; idx<=max_steps; idx++)   // max_steps = f(dt_mech, max_time)
        {
            update_motility_vector( dt_mech );
            // velocity += motility_vector;
            velocity[0] += motility_vector[0];
            velocity[1] += motility_vector[1];
            // std::cout <<"motility vector= "<<cell_xpos[idx]<<","<<cell_ypos[idx]<<std::endl;
            update_position( dt_mech );
            // std::cout <<idx<<") cell pos= "<<position[0]<<","<<position[1]<<std::endl;
            // if ((idx % 1000) == 0)  // equiv of <interval units="min">100</interval>
            if ((idx % idx_output) == 0)  // equiv of <interval units="min">100</interval>
            {
                // printf("%d,%d,%.4f,%.4f\n", id_run, run_num, position[0],position[1]);
                // std::cout <<idx/10<<","<< run_num<<","<< position[0]<<","<<position[1] << std::endl;
                // std::cout <<idx<<","<< run_num<<","<< position[0]<<","<<position[1] << std::endl;
                std::cout <<idx<<","<< run_num<<","<< position[0]<<","<<position[1] << ",222.34,52.86" <<std::endl;
            }
        }
        // printf("The number is %d and Pi is %.2f\n", num, pi);
        // printf("%d,%d,%.4f,%.4f\n", id_run, run_num, position[0],position[1]);
    }

    // std::cout <<"\n-------------"<<std::endl;
    // for (int idx=0; idx<cell_xpos.size(); idx++)
    // {
    //     std::cout <<"motility vector= "<<cell_xpos[idx]<<","<<cell_ypos[idx]<<std::endl;
    //     // std::cout <<"cell pos= "<<position[idx]<<","<<position[idx]<<std::endl;
    // }
}
