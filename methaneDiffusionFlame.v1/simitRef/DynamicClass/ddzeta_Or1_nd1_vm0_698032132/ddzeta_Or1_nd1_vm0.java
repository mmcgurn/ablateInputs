import java.io.Serializable;
import simit.numerics.interpolant.*;
public class ddzeta_Or1_nd1_vm0 implements PolyInterpolantBSdynamicBuildInterface, Serializable{
public double ddzeta(double[][][] phi, double[] zeta, double[][] zetagrid, int[] ijklow, int nd, double[][] _phi_interp){
double result = 0;
	 _phi_interp[0][0] = phi[0][0][0+ijklow[0]];
	 _phi_interp[0][1] = phi[0][0][1+ijklow[0]];
{
result= 0.;
double weight = 0.;
		weight = 0.;
		{
		 double prod = 1.;
		 prod = 1.;
		 prod = prod /(zetagrid[0][0] - zetagrid[0][1]); 
		 weight += prod;
		}
	result+= weight * _phi_interp[0][0];
		weight = 0.;
		{
		 double prod = 1.;
		 prod = 1.;
		 prod = prod /(zetagrid[0][1] - zetagrid[0][0]); 
		 weight += prod;
		}
	result+= weight * _phi_interp[0][1];
}
return result;
}
}
