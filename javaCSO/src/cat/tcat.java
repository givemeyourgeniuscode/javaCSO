package cat;

public class tcat {
public	double[] d_pos;			  // Record the cat's coordinate.
public	double val;				    // Fitness value, which is evaluated by the benchmark function.
public	int STFlag;				  // The flag controls the cat to move into the seeking mode or the tracing mode.
public	double[] d_vel;			  // The velocity corresponding to the position.
public tcat(int i)
{
	d_pos = new double[i];
    d_vel = new double[i];
}

}
