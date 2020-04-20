package CSO;

import cat.tcat;

import java.io.*;
import java.util.ArrayList;

import Benchmark.Benchmark;
public class CSO {

	// Parameter setting:

	/* Parameters for the Seeking Mode - Start */
	static int    SMP =           5;		    // The maximum contineuously stay and view times.
	static double SRD =		    0.2;     // The persentage of the perturbation on the selected dimensions in the "Tracing" mode.
	static double CDC =			0.8;     // Mutation Range form this point to the other.
	static int    SPC =	    	  1;       // Self-position Considering in the "Seeking" mode. Set "1" to activate; set "0" to deactivate.
	/* Parameters for the Seeking Mode - End   */

	/* Parameter for the Tracing Mode - Start */
	static double C1 =      2.0;     // Parameter for the velocity
	/* Parameter for the Tracing Mode - End   */

	/* Parameter for repicking the cats to move by the tracing mode in the next iteration. - Start */
	static double  MR =			0.02;	  // The ratio of cats in the "Seeking" mode.
	/* Parameter for repicking the cats to move by the tracing mode in the next iteration. - End   */

	static int i_Iteration =    60;
	
	int i_PSize;
	int i_Dim;
	int i_FN;
	double d_InitL;
	double d_InitR;
	double d_MaxVel;
	boolean b_Maximize;
	ArrayList<tcat> t_cat = new ArrayList<tcat>();
	
	CSO( int i_FN,int i_PSize, int i_Dim, double d_InitL, double d_InitR, double d_MaxVel, boolean b_Maximize )
	{
		this.i_FN=i_FN;
		this.i_PSize=i_PSize;
		this.i_Dim=i_Dim;
		this.d_InitL=d_InitL;
		this.d_InitR=d_InitR;
		this.d_MaxVel=d_MaxVel;
		this.b_Maximize=b_Maximize;
		
		for (int k = 0; k < (i_PSize + 1); k++)					     
		{
			
			 t_cat.add( new tcat(i_Dim));
		}
	}
	
	
public static void main (String[] args )
{
    CSO acso = new CSO(4 , 1 , 2 , -3 , 6 , 4 , true);
	int cycle = 5;
	long average=0;
	long[] timer = new long[cycle];
	for(int i=0;i<cycle;i++)
	{    
		//start timing
		long start =  System.currentTimeMillis();
		
		acso.CSO_Run();
		//end timing
		long end =  System.currentTimeMillis();
		timer[i]=(end-start)/1000;
	}
	for(int i=0;i<cycle;i++)
		average+=timer[i];
	average/=50;
	System.out.println("50次平均用时:  "+average);
}
	
	
	
	
double fabs(double n)
	{
		if (n>=0 )
			return n;
		else return -n;
	}
	
double Rand() // This subroutine returns a random value in the range [0.0, 1.0].
	{
		double d_ret;
		d_ret = ((double)Math.random() / (Math.pow(2.0, 31.0) - 1.0));
		return d_ret;
	}
	
void RandSeq(int i_num, int i_seq[]) // This subroutine returns the random sequence with i_num elements.
	{
		int[] i_tmp=new int[2];
		int i;
		double d_tmp;

		for (i = 0; i < i_num; i++)
			i_seq[i] = i;

		for (i = 0; i < i_num; i++)
		{
			d_tmp = (Rand()*(double)i_num);
			i_tmp[0] = (int)d_tmp;
			i_tmp[1] = i_seq[i_tmp[0]];
			i_seq[i_tmp[0]] = i_seq[i];
			i_seq[i] = i_tmp[1];
		}
	}


void CSO_Initialization( )
{
	int i, j, i_tmp;
	int[] i_seq;
	double d_tmp;
	double[] d_copy;

	int i_NTracing;                       // Contains the number of cats to be selected to move in the tracing mode.

	d_tmp = MR * (double)i_PSize;             // MR is defined in the "CSO_lib.h".
	i_NTracing = (int)d_tmp;                // The number of cats to move in the tracing mode.

	d_copy = new double[i_Dim];

	i_seq = new int[i_PSize];

	/* Randomly spread the cats into the solution space - Start */
	for (i = 0; i < (i_PSize + 1); i++)            // The for-loop for initializing the cats.
	{
		t_cat.get(i).STFlag = 0;                  // Set the Seeking/Tracing flag to be deactivated.
		for (j = 0; j < i_Dim; j++)                // The for-loop for initializin the dimensions and the velocities of a cat.
		{
			t_cat.get(i).d_pos[j] = ((Rand()*(d_InitR - d_InitL)) + d_InitL); // Set the coordinate of the cat randomly in the initial range.
			t_cat.get(i).d_vel[j] = ((Rand()*(2.0*d_MaxVel)) - d_MaxVel);   // Set the velocity of the cat randomly in the range of the maximum velocity.
		}                                   // End of the for-loop of the dimension setting.
	}                                     // End of the for-loop of the initialization.
	/* Randomly spread the cats into the solution space - End   */

	/* Calculate the randomly generated near best solution's fitness vaule - Start */
	for (j = 0; j < i_Dim; j++)
		d_copy[j] = t_cat.get(i_PSize).d_pos[j];  // Copy the coordinate of the cat for further usage.

	t_cat.get(i_PSize).val = Benchmark.Fuction(i_FN, i_Dim, d_copy);  // Calculate the fitness value of the randomly generated near best solution.
	/* Calculate the randomly generated near best solution's fitness vaule - End   */

	/* Pick i_NTracing cats and set their STFlage to be activated - Start */
	RandSeq(i_PSize, i_seq);

	for (j = 0; j < i_NTracing; j++)
		t_cat.get(i_seq[j]).STFlag = 1;
	/* Pick i_NTracing cats and set their STFlage to be activated - End   */

	
}

	
void CSO_Tracing_Mode(int i_tg) throws IOException
	{
		int i;

		double w;
		//ofstream fout("MatrixW.txt");
		//OutputStream fout = new FileOutputStream("MatrixW.txt");
		double d_minus = (0 - 1.0*d_MaxVel);

		for (i = 0; i < i_Dim; i++)
		{
			w = Rand();
			w *= 1000;
			// fout << w << endl;
          // fout.write((int) w);
          
			t_cat.get(i_tg).d_vel[i] = w * t_cat.get(i_tg).d_vel[i] + Rand()*C1*(t_cat.get(i_PSize).d_pos[i] - t_cat.get(i_tg).d_pos[i]); // Update the velocity of the i_tg cat.

			/* Check if the new velocity is larger/smaller than the limitation - Start */
			if (t_cat.get(i_tg).d_vel[i] > d_MaxVel)
				t_cat.get(i_tg).d_vel[i] = d_MaxVel;
			else if (t_cat.get(i_tg).d_vel[i] < d_minus)
				t_cat.get(i_tg).d_vel[i] = d_minus;
			/* Check if the new velocity is larger/smaller than the limitation - End   */

			t_cat.get(i_tg).d_pos[i] += t_cat.get(i_tg).d_vel[i]; // Update the coordinate of the i_tg cat.
		}
		// fout.close();
	}
	

void CSO_Evaluation()
{
	int i, j;
	double[] d_copy;

	d_copy = new double[i_Dim];

	for (i = 0; i < i_PSize; i++)
	{
		for (j = 0; j < i_Dim; j++)
		{
			d_copy[j] = t_cat.get(i).d_pos[j];
		}
		t_cat.get(i).val = Benchmark.Fuction(i_FN, i_Dim, d_copy);  // Calculate the fitness value of the ith cat.

		/* If the current solution is better than the stored one, update the best solution in the memory. - Start */
		if (b_Maximize)  // If the goal is to find the near maximum solution, goes into this part of the program.
		{
			if (t_cat.get(i_PSize).val < t_cat.get(i).val)
			{
				t_cat.get(i_PSize).val = t_cat.get(i).val;
				for (j = 0; j < i_Dim; j++)
				{
					t_cat.get(i_PSize).d_pos[j] = t_cat.get(i).d_pos[j];
				}
				//        printf("The near best solution is updated!\n");
			}
		}
		else            // If the goal is to find the near minimum solution, goes into this part of the program.
		{
			if (t_cat.get(i_PSize).val > t_cat.get(i).val)
			{
				t_cat.get(i_PSize).val = t_cat.get(i).val;
				for (j = 0; j < i_Dim; j++)
				{
					t_cat.get(i_PSize).d_pos[j] = t_cat.get(i).d_pos[j];
				}
				//        printf("The near best solution is updated!\n");
			}
		}
		/* If the current solution is better than the stored one, update the best solution in the memory. - End   */
	}
	
}
	

void CSO_Seeking_Mode(int i_tg)
{
	int i, j, i_tmp, i_SK, i_ModDm;
	double d_tmp, d_max, d_min;
	double[] d_Seeking ;             // The variable to contain the copy of the seeking
	double[] d_copy ,d_cpfns ,d_stfns ;
	int[] i_seq ;

	d_Seeking = new double[(SMP*i_Dim)];// Allocate the memory to generate SMP copies of the current coordinate.
	d_copy = new double[i_Dim];          // Allocate the memory to contain the coordinate, which will be sent to the benchmark function.
	d_cpfns = new double[SMP];           // Allocate the memory to contain the new fitness value.
	d_stfns = new double[SMP];           // Allocate the memory to contain the copied fitness values for sorting.
	i_seq = new int[i_Dim];              // Allocate the memory to contain the random sequence.

	d_tmp = CDC * (double)i_Dim;
	i_ModDm = (int)d_tmp;                   // The number of dimensions to be changed for one candidate.

	/* Make SMP copies of the current coordinate - Start */
	for (i = 0; i < SMP; i++)
		for (j = 0; j < i_Dim; j++)
			d_Seeking[((i*i_Dim) + j)] = t_cat.get(i_tg).d_pos[j];
	/* Make SMP copies of the current coordinate - End   */

	if (SPC==1)                                         // The current coordinate is also considered as a candidate.
	{
		i_SK = (SMP - 1);
		d_cpfns[(SMP - 1)] = t_cat.get(i_tg).val;
		d_stfns[(SMP - 1)] = d_cpfns[(SMP - 1)];
	}
	else                                            // The current coordinate is not considered as a candidate.
		i_SK = SMP;

	/* Generate the seeking spots as the candidates - Start */
	for (i = 0; i < i_SK; i++)                             // The for-loop for creating the seeking spots as the candidates.
	{
		RandSeq(i_Dim, i_seq);
		for (j = 0; j < i_ModDm; j++)
		{
			d_tmp = d_Seeking[((i*i_Dim) + i_seq[j])] * SRD;  // Modify the value on the selected dimension with SRD percent.
			if (Rand() < 0.5)                              // The change that the change becomes larger or smaller is equal.
				d_tmp *= (-1.0);
			d_Seeking[((i*i_Dim) + i_seq[j])] += d_tmp;     // Modify the value on the selected dimension.
		}
		for (j = 0; j < i_Dim; j++)
			d_copy[j] = d_Seeking[((i*i_Dim) + j)];         // Copy the coordinate of the seeking spot. 
		d_cpfns[i] = Benchmark.Fuction(i_FN, i_Dim, d_copy);      // Calculate the fitness value of the seeking spot.
		d_stfns[i] = d_cpfns[i];                        // Copy the fitness value for sorting.
	}                                               // End of the for-loop for creating the seeking spots.
	/* Generate the seeking spots as the candidates - End   */

	/* Calculate the chance to movie to the candidate spots - Start */
	for (i = 0; i < SMP; i++)                              // Bubble sort to sort the fitness value of the candidates - Start
	{
		for (j = 1; j < SMP; j++)
		{
			if (d_stfns[j] < d_stfns[j - 1])                 // The smaller element will be moved to the front of the sequence.
			{
				d_tmp = d_stfns[j - 1];
				d_stfns[j - 1] = d_stfns[j];
				d_stfns[j] = d_tmp;
			}
		}
	}                                               // Bubble sort to sort the fitness value of the candidates - End
	d_max = d_stfns[(SMP - 1)];                         // Record the maximum fitness value of the candidates.
	d_min = d_stfns[0];                               // Record the minimum fitness value of the candidates.

	for (i = 0; i < SMP; i++)                              // Calculate the chance of the candidate to be picked.
	{
		if ((d_stfns[0] - d_stfns[(SMP - 1)]) == 0.0)        // Check if all the candidates present the same fitness vaule.
			d_stfns[i] = 1.0;                             // If the candidates persent the same fitness, set the chance to be "1.0" for all of them since there is no difference to pick any of them.
		else
		{
			if (b_Maximize)                              // If the goal is to find the near maximum solution, FSb=FSmin, otherwise FSb=FSmax.
				d_stfns[i] = (fabs(d_cpfns[i] - d_min) / (d_max - d_min));
			else
				d_stfns[i] = (fabs(d_cpfns[i] - d_max) / (d_max - d_min));
		}
	}
	/* Calculate the chance to movie to the candidate spots - End   */

	/* Pick a candidate spot to move to - Start */
	d_tmp = Rand();                                   // Generate a random number in [0,1].
	for (i = 0; i < SMP; i++)
	{
		if (b_Maximize)
		{
			if (d_tmp <= d_stfns[i])
			{
				for (j = 0; j < i_Dim; j++)
					t_cat.get(i_tg).d_pos[j] = d_Seeking[(i*i_Dim) + j];
				break;
			}
		}
		else
		{
			if (d_tmp <= d_stfns[i])
			{
				for (j = 0; j < i_Dim; j++)
					t_cat.get(i_tg).d_pos[j] = d_Seeking[(i*i_Dim) + j];
				break;
			}
		}
	}
	/* Pick a candidate spot to move to - End   */

	
}	
	

void CSO_Movement( )
{
	int i_NTracing;                       // Contains the number of cats to be selected to move in the tracing mode.
	double d_tmp;
	int i;
	int[] i_seq;
	d_tmp = MR * (double)i_PSize;             // MR is defined in the "CSO_lib.h".
	i_NTracing = (int)d_tmp;                // The number of cats to move in the tracing mode.

	i_seq = new int[i_PSize];             // Allocate a memory space to contain the series number of the cats.

	for (i = 0; i < i_PSize; i++)                // All the cats should take turns to move according to their STFlags.
	{
		if (t_cat.get(i).STFlag==1)
			try {
				CSO_Tracing_Mode(i);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				System.out.println("文件读写有误");
				e.printStackTrace();
			}
		else                                // This cat should go into the Seeking mode.
			CSO_Seeking_Mode( i);

		t_cat.get(i).STFlag = 0;                  // Set the STFlags of all the cats to be deactivated.
	}

	/* Pick i_NTracing cats and set their STFlage to be activated - Start */
	RandSeq(i_PSize, i_seq);

	for (i = 0; i < i_NTracing; i++)
		t_cat.get(i_seq[i]).STFlag = 1;
	/* Pick i_NTracing cats and set their STFlage to be activated - End    */

}


void CSO_Run()
{
   						          // The populations are called "t_cat" in this program.
	int i, j, k;
					    
	

	CSO_Initialization();
		
	for (j = 0; j < i_Iteration; j++)				      // The for-loop for repeating the iteration.
		{
			CSO_Evaluation(); // Evaluate the fitness value of every cat in the current iteration.
			CSO_Movement();  // Move the cats in the solution space.
			
			if ((j == 0) || (((j + 1) % (i_Iteration/4)) == 0))
				System.out.println("[Iteration]: ["+(j + 1)+ "]=" +t_cat.get(i_PSize).val);

		}										                    // End of the for-loop of Iteration.

	
}
	
	
	
}
