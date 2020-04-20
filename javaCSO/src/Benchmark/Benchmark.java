package Benchmark;

public class Benchmark {
private static double d_PIE=3.141592653589793238462643383279;	// Defines the ratio of the circumference.
private static double d_PII=6.283185307179586476925286768758;	// Defines twice of the ration of the circumference.
private static int  i_FuncNum=3;					// i_FuncNum defines how many test functions are built in the Benchamrk library.
Benchmark() 
{
	
}

public static double  Fuction(int i_fn, int i_dim, double d_pos[])
{
	int i, j;
	double[] d_tmp=new double[4];	// variables for temporary storage.
	double d_result = 0;	// cotains the result to return to the main program.

	/* initialize the parameters - Start */
	for (i = 0; i < 4; i++)
		d_tmp[i] = 0.0;
	/* initialize the parameters - End   */

	if (i_fn == 1)	// Rosenbrock Function
	{
		for (j = 1; j < (i_dim - 1); j++)
			d_tmp[0] += (100.0*(Math.pow((d_pos[j] - Math.pow(d_pos[j - 1], 2.0)), 2.0)) + Math.pow((d_pos[j - 1] - 1.0), 2.0));

		d_result = d_tmp[0];
	}
	else if (i_fn == 2)	// Rastrigrin Function
	{
		for (j = 0; j < i_dim; j++)
		{
			d_tmp[0] = Math.pow(d_pos[j], 2.0);
			d_tmp[1] = (10.0*Math.cos(d_PII*d_pos[j]));
			d_tmp[2] += (d_tmp[0] - d_tmp[1] + 10.0);
		}
		d_result = d_tmp[2];
	}
	else if (i_fn == 3)	// Griewank Function
	{
		d_tmp[2] = 1.0;

		for (j = 0; j < i_dim; j++)
		{
			d_tmp[0] += (Math.pow((d_pos[j] - 100.0), 2.0));

			d_tmp[1] = ((d_pos[j] - 100.0) / Math.sqrt((j + 1)));
			d_tmp[2] *= Math.cos(d_tmp[1]);
		}
		d_result = (0.00025*d_tmp[0] - d_tmp[2] + 1.0);
	}
	
		else if(i_fn==4)	// Test Function 4
		{
			
			return 25-Math.pow(d_pos[0]-1,2)-Math.pow(d_pos[1]-2, 2);
			//-4<x<6,-3<y<7
		}/*
		else if(i_fn==5)  // Test Function 5
		{
		}
		else // Test Function 6
		{
		}
	*/

	return d_result;
}

}
