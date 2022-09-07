package simit.domain.hydro_domain.hydro_region.process.combustion.manifold_generation.flamelet_generator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import simit.domain.hydro_domain.hydro_material.Compressible_Material;
import simit.domain.hydro_domain.hydro_material.Compressible_Material_Soot;
import simit.domain.hydro_domain.hydro_mesh.BS_Mesh;
import simit.domain.hydro_domain.hydro_mesh.Cell_Centered_BS_Mesh;
import simit.domain.hydro_domain.hydro_mesh.OneDSphericallySymmetricMesh;
import simit.domain.hydro_domain.hydro_mesh.Vertex_Centered_BS_Mesh;
import simit.domain.hydro_domain.hydro_region.process.Process;
import simit.domain.hydro_domain.hydro_region.process.chemical_rxn.Chemical_Equilibrium;
import simit.domain.hydro_domain.hydro_region.process.chemical_rxn.Chemical_RXN;
import simit.domain.hydro_domain.hydro_region.process.combustion.manifold_generation.RPV_CxHyOz_Flamelet;
import simit.domain.hydro_domain.hydro_region.process.molecular.oneD_Cartesian.Diffusion_mflux_oneD_Cartesian_Master;
import simit.domain.hydro_domain.hydro_region.process.molecular.oneD_Cartesian.sidesets.oneDsolid_Yeq_Diffusion_SideSet;
import simit.domain.hydro_domain.hydro_region.process.molecular.twoD_Cartesian.sidesets.Fixed_Dirichlet;
import simit.domain.hydro_domain.hydro_region.sidesets.Hydro_Region_SideSet;
import simit.domain.hydro_domain.hydro_region.sidesets.Hydro_Region_SideSet_Manager;
import simit.domain.particle_domain.liquids.Methanol_Fuel;
import simit.domain.particle_domain.liquids.Wax_Fuel;
import simit.domain.particle_domain.liquids.nHeptane_Fuel;
import simit.domain.particle_domain.liquids.nHexadecane_Fuel;
import simit.domain.particle_domain.particle_types.Hydrocarbon.HCM.Hydrocarbon_Model;
import simit.domain.particle_domain.particle_types.Hydrocarbon.HCM.Multi_zone;
import simit.domain.particle_domain.particle_types.Hydrocarbon.HCM.Two_zone;
import simit.domain.radiation_domain.RayTracing.EmpericalGasEmissivityModels.GasAbsorptionModel;
import simit.domain.radiation_domain.RayTracing.EmpericalGasEmissivityModels.ZimmerGasModel;
import simit.domain.radiation_domain.RayTracing.MeanSootModels.SootMeanGasAbsorption;
import simit.domain.radiation_domain.RayTracing.MeanSootModels.SootRosselandMeanAbsorption;
import simit.domain.radiation_domain.RayTracing.Process.RayTracingFull1DPlanar;
import simit.io.exceptions.NumericalException;
import simit.io.exceptions.ParticleException;
import simit.io.exceptions.SMTException;
import simit.io.rtesupport.RTE;
import simit.io.rtesupport.Runner;
import simit.numerics.Function;
import simit.numerics.constants.Const;
import simit.numerics.geometry.SphericalDistribution;
import simit.numerics.primatives.FLAGS;
import simit.numerics.primatives.MassFrac;
import simit.numerics.primatives.Pressure;
import simit.numerics.primatives.Property;
import simit.numerics.primatives.Temperature;
import simit.numerics.primatives.Velocity;
import simit.numerics.primatives.vec;
import simit.numerics.tables.GeneralTable;
import simit.shared_resources.Material;
import simit.shared_resources.SMTarray;
import simit.shared_resources.SMTstring;
import simit.shared_resources.var;
import simit.shared_resources.eos.EOS_Table;
import simit.shared_resources.eos.Ideal_EOS;
import simit.shared_resources.eos.Ideal_EOS_LHF;
import simit.shared_resources.eos.Ideal_Gas_EOS;
import simit.shared_resources.eos.spec;
import simit.shared_resources.rxns.Kinetics_Model;
import simit.shared_resources.rxns.Kinetics_Model_Parser;
import simit.shared_resources.rxns.YamlParser;
import simit.shared_resources.transport_prop.KineticTheoryTransport;
import simit.shared_resources.transport_prop.SutherlandTransport;
import simit.shared_resources.transport_prop.TransportPropsModel;

/**
 * Purpose of this class is to solve steady-state flames starting with a relatively large domain defined by Lref. The domain is decreased by a factor
 * scalefac(=0.5) for every new flame. For every set of BCs, a family of flames are saved to an xml formatted file. The BCs are read in as xml formatted files.
 * The number of flames per set of BCs is defined in the main statement.
 * 
 * @author ped3
 *
 */
public class oneDflame_SS_Solver_Master implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	oneDflame_Region region;
	final static int NX = 100;
	final static double dto = 1.e-7;
	final static double scalefac = 0.9; // scale factor to scale each new flame based on previous
	final double Lref;
	final oneDflames flames;
	Cell_Centered_BS_Mesh cmesh;
	Map<String, Property> BC1, BC2;
	Chemical_RXN rxn;
	private boolean SootToggle;
	private static boolean RadiationToggle = false;
	private static boolean CoupledBCToggle = true;
	private static boolean SphericalToggle = false;
	private static double Sphericalr0 = 0; // Initial droplet radius
	private static double SphericalScale = 12; // Initial droplet domain Scaling

	// Store the type of Initializer
	public static enum Initializer {
		Equillibrium, Linear, Flamelet, Multizone
	};

	public static Initializer InitMethod = Initializer.Equillibrium;

	private static oneDflames InitializingFGM = null;
	private static int numThreads = 1;
	private static int dtChem = 1;

	public oneDflame_SS_Solver_Master(String name, RTE rte, String mech_name, Map<String, Property> BC1, Map<String, Property> BC2, double Lref) throws Exception {
		this(name, rte, mech_name, BC1, BC2, Lref, null);
	}

	public oneDflame_SS_Solver_Master(String name, RTE rte, String mech_name, Map<String, Property> BC1, Map<String, Property> BC2, double Lref, oneDflames InitializingFGM) throws Exception {
		oneDflame_SS_Solver_Master.InitializingFGM = InitializingFGM;
		FLAGS.VERBOSE_IT = true;
		this.BC1 = BC1;
		this.BC2 = BC2;
		this.Lref = Lref;
		Kinetics_Model kinetics;
		String[] species;
		Ideal_EOS eosOg;
		Ideal_EOS eos;
		Ideal_EOS eosTable;
		Ideal_EOS.Pref_for_ConstP = Const.atm;
		// Load in the mech
		Chemical_Equilibrium CE = new Chemical_Equilibrium(mech_name);
		if (mech_name.contains(".yaml")) {
			YamlParser parser = new YamlParser(mech_name);
			kinetics = parser.getKinetics();
			species = kinetics.get_Species();
			eosOg = parser.getEOS();
			eosOg.CONSTPTOG = true;
			eosTable = new EOS_Table(eosOg, SMTstring.addOnlyNewStrings(species, CE.species));
			eosTable.CONSTPTOG = true;
		} else {
			kinetics = new Kinetics_Model_Parser(mech_name);
			species = kinetics.get_Species();
			// Add any other species needed from the mech
			String[] allSpecies = SMTstring.addOnlyNewStrings(species, CE.species);
			eosOg = new Ideal_Gas_EOS(allSpecies, mech_name);
			eosOg.CONSTPTOG = true;
			eosTable = new EOS_Table(eosOg, allSpecies);
			eosTable.CONSTPTOG = true;
		}

		this.CheckForSoot(species);
		if (SootToggle)
			eos = new Ideal_EOS_LHF(species, eosTable, new String[] { "C(s)" }, new double[] { Const.rhoC });
		else
			eos = eosTable;
		eos.CONSTPTOG = true;
		((MassFrac) BC1.get(var.Yi)).add(species, 0.); // expand BCs to include all species in reaction set
		((MassFrac) BC2.get(var.Yi)).add(species, 0.); // if it doesn't already contain them.

		TransportPropsModel transport = new SutherlandTransport();

		if (SphericalToggle) {
			cmesh = new Cell_Centered_BS_Mesh(new Vertex_Centered_BS_Mesh(new int[] { NX }, new vec(Sphericalr0), new vec(Sphericalr0 * SphericalScale)));
			cmesh = new OneDSphericallySymmetricMesh(cmesh, new SphericalDistribution(1, 1));
		} else
			cmesh = new Cell_Centered_BS_Mesh(new Vertex_Centered_BS_Mesh(new int[] { NX }, new vec(0.), new vec(Lref)));

		Compressible_Material gas;
		if (SootToggle)
			gas = new Compressible_Material_Soot(eos, transport, cmesh, species);
		else
			gas = new Compressible_Material(eos, transport, cmesh, species);

		gas = this.InitializeGas(gas, CE, species, eos, transport);

		Vector<Process> proc = new Vector<Process>();

		// Create a chemical rxn spec
		Chemical_RXN.ChemicalRnxSpec spec = new Chemical_RXN.ChemicalRnxSpec();
		spec.setMassDtFrac(dtChem);
		spec.setEnerDtFrac(dtChem);
		spec.setRTol(1.E-8);
		spec.setATol(1.E-8 * 1.E-2);
		spec.setIntegrationConstraint(Chemical_RXN.IntegrationConstraint.MassEner);
		spec.setWithJVODE(true);
		spec.setNumberOfThreads(numThreads);
		spec.setOneAtmPressure(true);
		this.rxn = new Chemical_RXN(gas, kinetics, spec);
		rxn.output();
		proc.add(rxn);

		// Add Radiation `Solver?
		RayTracingFull1DPlanar RT1D = null;
		if (RadiationToggle) {
			GasAbsorptionModel gasabs = new ZimmerGasModel(species);
			SootMeanGasAbsorption sootAbs = new SootRosselandMeanAbsorption(7.);
			RT1D = new RayTracingFull1DPlanar(gas, "RT1DSolver", gasabs, sootAbs, 90, new double[] { 1., 1. });
		}
		if (RT1D != null)
			proc.add(RT1D);

		Diffusion_mflux_oneD_Cartesian_Master.implicit_tog = true;
		Hydro_Region_SideSet bc1, bc2;
		if (!SphericalToggle) {
			bc1 = new Fixed_Dirichlet(0, gas, BC1);
			if (!CoupledBCToggle)
				bc2 = new Fixed_Dirichlet(1, gas, BC2);
			else
				bc2 = new oneDsolid_Yeq_Diffusion_SideSet(1, gas, BC2, RT1D);
		} else {
			if (CoupledBCToggle)
				bc1 = new oneDsolid_Yeq_Diffusion_SideSet(0, gas, BC2, RT1D);
			else
				bc1 = new Fixed_Dirichlet(0, gas, BC2);
			bc2 = new Fixed_Dirichlet(1, gas, BC1);
		}

		Hydro_Region_SideSet_Manager DSSM = new Hydro_Region_SideSet_Manager(new Hydro_Region_SideSet[] { bc1, bc2 }, gas.mesh);
		proc.add(new Diffusion_mflux_oneD_Cartesian_Master(gas, DSSM, SphericalToggle));
		//// Create Mixture Fraction Calculator
		MixFrac_Calc MFCalc = new MixFrac_Calc(((MassFrac) BC2.get(var.Yi)).getValues(), ((MassFrac) BC1.get(var.Yi)).getValues(), species, mech_name);
		this.region = new oneDflame_Region(proc, gas, MFCalc);

		flames = new oneDflames(name, rxn, gas, new File(RTE.dir.getResultDir(), "oneDflames"), "oneDflames");
	}

	private Compressible_Material InitializeGas(Compressible_Material gas, Chemical_Equilibrium CE, String[] species, Ideal_EOS eos, TransportPropsModel transport) throws Exception {
		System.out.println("=================================");
		System.out.println("   Initializing Original Guess   ");
		System.out.println("=================================");
		switch (oneDflame_SS_Solver_Master.InitMethod) {
		case Equillibrium: {
			// BC's flopped in spherical domain
			if (!SphericalToggle)
				Material.init_region(new Temperature(new linearfunc(cmesh, BC1.get(var.T), BC2.get(var.T))), new Pressure(Const.atm), new Velocity(0.), new MassFrac(species, new linearfunc(cmesh, BC1.get(var.Yi), BC2.get(var.Yi))), gas);
			else
				Material.init_region(new Temperature(new linearfunc(cmesh, BC2.get(var.T), BC1.get(var.T))), new Pressure(Const.atm), new Velocity(0.), new MassFrac(species, new linearfunc(cmesh, BC2.get(var.Yi), BC1.get(var.Yi))), gas);
			for (int i = cmesh.isc; i <= cmesh.iec; i++) {
				CE.solveEquil(gas.getYi(i, 0, 0), gas.T[0][0][i], gas.p[0][0][i], "HP");
				double[] Yigas = gas.eos.XitoYi(CE.Xi, CE.specID);
				for (int ns = 0; ns < gas.nspec; ns++)
					if (Yigas[ns] < Const.SMALL)
						Yigas[ns] = 0.; // clip anything very small to zero.
				MassFrac.normalize(Yigas);
				gas.set_point_TPvelYi(CE.T_out, Ideal_EOS.Pref_for_ConstP, gas.get_vel_vec(0, 0, i).getDoubleArray(), Yigas, null, 0, 0, i);
			}
			return gas;
		}
		case Linear: {
			double[] Yf = ((MassFrac) BC2.get(var.Yi)).getValues();
			double[] Yox = ((MassFrac) BC1.get(var.Yi)).getValues();
			int nO2 = spec.lookupID("O2");
			int nN2 = spec.lookupID("N2");
			double TLeft = ((Temperature) BC1.get(var.T)).getValue();
			double TRight = ((Temperature) BC2.get(var.T)).getValue();
			double[] Ypt = new double[species.length];
			for (int i = cmesh.is; i <= cmesh.ie; i++) {
				SMTarray.zero(Ypt);
				double TMax = 2200;
				double Peak = 2. / 10.;
				double x = gas.mesh.get_loc(0, 0, i).value(0);
				Ypt[nO2] = Yox[nO2] * (1 - x / Lref);
				double Sum = 0;
				for (int ns = 0; ns < species.length; ns++) {
					if (ns != nN2)
						Ypt[ns] += Yf[ns] * (x / Lref);
					Sum += Ypt[ns];
				}
				Ypt[nN2] = 1. - Sum;
				double Temp = 0;
				if (x <= Peak * Lref)
					Temp = TLeft + (TMax - TLeft) * (x / (Peak * Lref));
				if (x >= Peak * Lref)
					Temp = TMax + (TRight - TMax) * (x - Peak * Lref) / (Lref - Peak * Lref);
				gas.set_point_TPvelYi(Temp, Ideal_EOS.Pref_for_ConstP, gas.get_vel_vec(0, 0, i).getDoubleArray(), Ypt, null, 0, 0, i);
			}
			return gas;
		}
		case Flamelet: {
			// Flamelet Initialization
			Map<String, Map<String, double[]>> flamesConstr = oneDflame_SS_Solver_Master.InitializingFGM.get_flames(); // grabs the flames
			Map<String, double[]> flameInterp = flamesConstr.get(flamesConstr.keySet().iterator().next()); // grabs the first key
			double[] Temp = flameInterp.get(var.T);// grabs temperature values
			double[] X = flameInterp.get(var.pos);
			// Create Tables
			GeneralTable TempTable = new GeneralTable(X, Temp);
			GeneralTable[] SpeciesTable = new GeneralTable[species.length];
			// Need to grab all Yi values that are known
			double[] Yi;
			for (int ns = 0; ns < species.length; ns++) {
				Yi = flameInterp.get(var.Yi + species[ns]);
				if (Yi == null) {
					Yi = new double[X.length];
				}
				SpeciesTable[ns] = new GeneralTable(X, Yi);
			}
			// Now interpolate values
			double T;
			Yi = new double[species.length];
			for (int i = cmesh.is; i <= cmesh.ie; i++) {
				T = TempTable.lookUp(cmesh.get_loc(0, 0, i).value(0));
				for (int ns = 0; ns < species.length; ns++)
					Yi[ns] = SpeciesTable[ns].lookUp(cmesh.get_loc(0, 0, i).value(0));
				gas.set_point_TPvelYi(T, Ideal_EOS.Pref_for_ConstP, gas.get_vel_vec(0, 0, i).getDoubleArray(), Yi, null, 0, 0, i);
			}
			return gas;
		}
		case Multizone: {

			// INPUT variables:
			double yio2 = 0.23; 
			double Tp = 370;
			double Dp = 2*cmesh.x[0][0][0][0];
			double Tinf = Const.Tref;
			double Pinf = Const.atm;
			double n2Do2=3.76; // For CxHyOz initialization
			
			nHeptane_Fuel fuel = new nHeptane_Fuel();	
//			Wax_Fuel fuel = new Wax_Fuel();	
//			nHexadecane_Fuel fuel = new nHexadecane_Fuel();

			
			
			String[] speclabel = new String[] { fuel.formula()[0], "O2", "H2O", "CO2", "N2", "CO", "H2", "OH", "O2P" };
			spec.init(speclabel);

			Ideal_Gas_EOS eos2 = new Ideal_Gas_EOS(speclabel);
			RPV_CxHyOz_Flamelet flame2 = new RPV_CxHyOz_Flamelet(fuel.formula()[0], eos2,n2Do2);
			TransportPropsModel transport2 = new KineticTheoryTransport(eos2);
			Hydrocarbon_Model GPM;
			GPM = new Multi_zone(fuel, eos2, transport2, flame2);

			Map<String, Double> gas_farfield = new LinkedHashMap<String, Double>();
			
			gas_farfield.put(var.T, Tinf);
			gas_farfield.put(var.p, Pinf);

			LinkedHashMap<String, Double> YiZ0 = flame2.YiZ0();
			for (String key : YiZ0.keySet()) {
				gas_farfield.put(var.Yi + key, YiZ0.get(key));
				gas_farfield.put(var.Yi + "N2", 1 - yio2);
				gas_farfield.put(var.Yi + "O2", yio2);
			}
			Map<String, Double> gas_sol = new TreeMap<String, Double>();
			Map<String, double[][][]> gas_init = new TreeMap<String, double[][][]>();
			

			int nzones = 1;
			double[][] data = null;
			try {
//					GPM.solve(Dp, Tpold, gas_farfield,gas_sol);
				data = GPM.solve1DSphereMulti(Dp, Tp, gas_farfield, gas_init);
				nzones = data[0].length;
			} catch (NumericalException e) {
				System.out.println("Multizone shit itself");
				System.exit(1);
			}

			String mech = "GRIMech30";
			Kinetics_Model kinetics = null;
			try {
				kinetics = new Kinetics_Model_Parser(mech);
			} catch (FileNotFoundException e) {
				System.out.println("Unable to create Kinetics_Model_Parser in Chemical_Equilibrium");
				e.printStackTrace();
			}

			String[] species2 = kinetics.get_Species();
			Chemical_Equilibrium CE2 = new Chemical_Equilibrium(mech);
			double[][] Yitemp1 = new double[species.length][nzones];
			double[] Yifuel = new double[species2.length];

			// FIX THIS

			for (int n = 0; n < CE2.specID.length; n++) {
				if (CE2.species[n].equals("C"))
					Yifuel[n] = fuel.m[0] * Const.MWC / fuel.MW();

				if (CE2.species[n].equals("H"))
					Yifuel[n] = fuel.n[0] * Const.MWH / fuel.MW();

				if (CE2.species[n].equals("O"))
					Yifuel[n] = fuel.x[0] * Const.MWO / fuel.MW();

			}


			cmesh = new Cell_Centered_BS_Mesh(new Vertex_Centered_BS_Mesh(new int[] { nzones }, new vec(data[0][0] * Dp * 0.5), new vec(data[0][nzones-1] * Dp * 0.5)));

//				double[] Yifueleq = new double[] {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//				0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//				0.,0.,0.,0.,0.55,0.,0.,0.,0.,0.,
//				0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//				0.,0.,0.,0.,0.,0.,0.,0.,0.45,0.,0.,0.,0.};//heptane GRI
				
				double[] Yifueleq = new double[CE.species.length];
				Yifueleq[96]=1.;

//			double[] Yifueleq = new double[] { 2.6605252679403324e-08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.03130099022620085, 0., 0., 0., 0., 0., 0.0, 0.0, 0., 0.019041133292868283, 0.0, 0.9085604652430157, 0., 0.012981434565768965, 0.,
//					0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0., 0.02811595000348999, 0., 0. };// wax
				
				
				//Setting up species
//			for (int i = 0; i < nzones; i++) {
//
////					double T=data[1][i];
////					double P=data[2][i];
////					Yifueleq=CE2.mapSol1D(Yifuel, T, P);
//				
//				String tempstring;
//				for (int n = 0; n < CE2.specID.length; n++) {
//					boolean a = true;
//					tempstring=CE.species[n];
//					for (int q = 0; q < speclabel.length; q++) {
//						if (CE.species[n].equals(speclabel[q])) {
//							Yitemp1[n][i] = Yifueleq[n] * data[3][i] + data[q + 3][i];
//							a = false;
//						}
//					}
//					if (a)
//						Yitemp1[n][i] = Yifueleq[n] * data[3][i];
//				}
//			}
			
			for (int i = 0; i < nzones; i++) {

//				double T=data[1][i];
//				double P=data[2][i];
//				Yifueleq=CE2.mapSol1D(Yifuel, T, P);
			
			String tempstring;
			for (int n = 0; n < CE2.specID.length; n++) {
				boolean a = true;
				tempstring=CE.species[n];
				for (int q = 0; q < speclabel.length; q++) {
					if (tempstring.toUpperCase().equals(speclabel[q])) {
						Yitemp1[n][i] = Yifueleq[n] * data[3][i] + data[q + 3][i];
						a = false;
					}
				}
				if (a)
					Yitemp1[n][i] = Yifueleq[n] * data[3][i];
			}
		}

			//mapping into vertex mesh
			
			
			
			for (int i = 0; i < nzones ; i++) {
				if (i==0)
					cmesh.vertex_mesh.x[0][0][0][i] = data[0][i] * Dp * 0.5 ;
				else if (i==(nzones-1))
					cmesh.vertex_mesh.x[0][0][0][i] = data[0][i] * Dp * 0.5 ;
				else 
					cmesh.vertex_mesh.x[0][0][0][i] = (data[0][i] + data[0][i + 1])/2 * Dp * 0.5  ;
			}
			cmesh = new OneDSphericallySymmetricMesh(cmesh, new SphericalDistribution(1, 1));
//				gas.updatemesh(cmesh);

			if (SootToggle)
				gas = new Compressible_Material_Soot(eos, transport, cmesh, species);
			else
				gas = new Compressible_Material(eos, transport, cmesh, species);
			double[] Temperature=zonal2mesh(data[1]);
			double[][] Yitemp2=new double[CE.specID.length][nzones+1];
			for (int n = 0; n < CE.specID.length; n++) {
				Yitemp2[n] =zonal2mesh(Yitemp1[n]);
			}
					
			double[] Yi = new double[CE.specID.length];
			for (int i = cmesh.is; i <= cmesh.ie; i++) {
				for (int n = 0; n < CE.specID.length; n++) {
					Yi[n] = Yitemp2[n][i];
				}				
				gas.set_point_TPvelYi(Temperature[i], Ideal_EOS.Pref_for_ConstP, gas.get_vel_vec(0, 0, i).getDoubleArray(), Yi, null, 0, 0, i);
			}
			return gas;
		}
		default: {
			System.out.println("Initializing method Unknown!");
			System.exit(1);
			return gas;

		}
		}

	}
	private double[] zonal2mesh(double[] Temp){

		double[] retmat=new double[Temp.length+1];
		for (int i = 0; i < Temp.length+1 ; i++) {
			if (i==0)
				retmat[i] = Temp[i];
			else if (i==(Temp.length))
				retmat[i] = Temp[i-1]  ;
			else 
				retmat[i] = (Temp[i] + Temp[i -1])/2;
		}
		return retmat;
	}
	private void CheckForSoot(String[] species) {
		this.SootToggle = SMTstring.contains(species, "C(s)");
	}

	class linearfunc implements Function {
		private static final long serialVersionUID = 1L;
		double xmax, xmin;
		double[] leftBC, rightBC;
		int size;

		public linearfunc(BS_Mesh mesh, Property leftBC, Property rightBC) {
			this.xmax = mesh.xmax[0];
			this.xmin = mesh.xmin[0];
			this.leftBC = leftBC.getValues();
			this.rightBC = rightBC.getValues();
			this.size = Math.min(leftBC.getValues().length, rightBC.getValues().length);
		}

		public linearfunc(BS_Mesh mesh, double leftBC, double rightBC) {
			this.xmax = mesh.xmax[0];
			this.xmin = mesh.xmin[0];
			this.leftBC = new double[] { leftBC };
			this.rightBC = new double[] { rightBC };
			this.size = 1;
		}

		@Override
		public double[] eval(double... pos) {
			double[] result = new double[size];
			for (int r = 0; r < size; r++) {
				result[r] = leftBC[r] + (rightBC[r] - leftBC[r]) * (pos[0] - xmin) / (xmax - xmin);
			}
			return result;
		}
	}

	public void solve(RTE rte, int nf) throws SMTException {
		Runner r = rte.runner();
		r.setExitState(Runner.ExitState.CONTINUE);
//		SaveStateManager.SERIALIZE_DATA = false;

		// cycle through flames - reducing the domain
		double L = Lref;
		for (int f = 0; f < nf; f++) {
			System.out.println("Solving for flame: " + f);
			System.out.println("----------------------");
			if (f == 0) {
				r.run(this.region, dto, true);
			} else { // new strained flame
				r.run(this.region, dto, false);
			}
			flames.add_flame(Double.toString(L));
			// re-scale the length
			L *= scalefac;
			Vertex_Centered_BS_Mesh vmesh = cmesh.vertex_mesh;
			for (int i = 0; i < vmesh.nx[0]; i++) {
				vmesh.x[0][0][0][i] = scalefac * vmesh.x[0][0][0][i];
			}
			cmesh.initVertexMesh(vmesh);
			cmesh.calcMetrics();
		}
		// finally add 1 more "flame" corresponding to clearly an extinguished flame assuming a linear distribution from the BCs.
		Material.init_region(new Temperature(new linearfunc(cmesh, BC1.get(var.T), BC2.get(var.T))), new Pressure(Const.atm), new Velocity(0.), new MassFrac(region.mat.speclabel, new linearfunc(cmesh, BC1.get(var.Yi), BC2.get(var.Yi))), region.mat);
		if (rxn.source == null)
			rxn.source = new double[this.region.nspeceq][this.cmesh.nx[2]][this.cmesh.nx[1]][this.cmesh.nx[0]]; // needed in limit when rxn processes isn't
																												// added - checking diffusion.
		SMTarray.zero(rxn.source);
		flames.add_flame(Double.toString(0.));
	}

	void save() {
		// finally write all flames to an XML formatted file.
		flames.write_xml();
	}

	public static void main(String[] args) throws Exception {
		// Uncomment If Running with a coupled Boundary or Radiation
		oneDflame_SS_Solver_Master.RadiationToggle = false;
		oneDflame_SS_Solver_Master.CoupledBCToggle = true;
		oneDflame_SS_Solver_Master.SphericalToggle = false;
		oneDflame_SS_Solver_Master.Sphericalr0 = .0005; // 1 mm Diameter
		oneDflame_SS_Solver_Master.SphericalScale = 12; // A little more than the flame stand off distance
//		Diffusion_mflux_oneD_Cartesian_Master.toggleThermoDiff = true;
		// Select Error Variable
		oneDflame_Region.ErrorVariableNames = new String[] { "T" };
		oneDflame_Region.NumErrorVariables = oneDflame_Region.ErrorVariableNames.length;
		oneDflame_Region.err_max = .01;
		// Set Dt Factors and Number of Threads
		oneDflame_SS_Solver_Master.dtChem = 30;
		Diffusion_mflux_oneD_Cartesian_Master.dtF_diffusion = 50;
		oneDflame_SS_Solver_Master.numThreads = 5;
		// outputs
		oneDflame_Region.IT_per_output = 1;
		oneDflame_Region.tstep_per_IT = 10;
//		oneDflame_Region.max_IT = 5;

		// Select The Initialization Method
		oneDflames ODFInit = null;
		oneDflame_SS_Solver_Master.InitMethod = oneDflame_SS_Solver_Master.Initializer.Linear;
//		ODFInit = new oneDflames(new File("C:\\Users\\klbud\\git\\SIMIT_WEB_FILES\\MANIFOLD\\MMA_Manifold"), "oneDflames_Twall=653.0");
		// Select Starting Length and number of flames
		double Lstart = 0.01;
		int nflames = 3; // solve for at least 3 flames - need this many to defined BS manifold.
		RTE rte = new RTE(args, "GRISphericalheptane_reduced" + File.separator + "oneDflame_steady_Lmax=" + Lstart);

		// Select Mechanism
		String mech_name = "GRIMech30";
//		String mech_name = "HeptaneReduced.yaml";
//		mech_name = "GRIMech30.yaml";
//		mech_name = "GRIMech30Soot";
//		mech_name = "GRIMech30Soot.yaml";
//		mech_name = "MMAReduced.yaml";
//		mech_name = "MMAReducedSoot.yaml";

//		String fuelBC = "fuelBCPMMA";
//		String folder = "HEPT_AIR_GRI_BC";
		String folder = "Wax_O2_BC";
		String fuelBC = "fuelBC";
		String oxBC = "oxBC";

		// Download All the required Files
		File f = RTE.file.remoteFile("MANIFOLD/" + folder, fuelBC + ".xml");
		f = RTE.file.remoteFile("MANIFOLD/" + folder, oxBC + ".xml");
		f = RTE.file.remoteFile("MANIFOLD", folder);

		// Read in the Boundary Conditions
		Map<String, Map<String, Property>> BCfuellist = Property.Reader(f, fuelBC);
		Map<String, Map<String, Property>> BCoxlist = Property.Reader(f, oxBC);

		// Get the start time
		long startTime = System.currentTimeMillis();

		for (String key : BCfuellist.keySet()) {
			oneDflame_SS_Solver_Master FS = new oneDflame_SS_Solver_Master(key, rte, mech_name, BCoxlist.get(key), BCfuellist.get(key), Lstart, ODFInit);
			FS.solve(rte, nflames);
			FS.save();
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.println("Total Time (min): " + elapsedTime * 0.001 / 60.);

		rte.shutdown();
	}
}
